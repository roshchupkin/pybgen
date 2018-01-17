import os
import struct
from collections import OrderedDict
import bitarray as ba
import zlib
import numpy as np
import pandas as pd
try:
    import zstd
except:
    print("No module ZSTD")


__author__ = "Gennady Roshchupkin"
__email__ = "gena.roshchupkin@gmail.com"
__license__ = "GNU General Public License v3"
__version__ = '0.1.0'

class Bgen_identifier(object):
    def __init__(self, bgen, N, Lh, offset):
        self.identifier_block_len = ba.bitarray(endian="little")
        self.identifier_block_len.fromfile(bgen, 4)
        self.identifier_block_len = struct.unpack("<L", self.identifier_block_len)[0]

        self.N = ba.bitarray(endian="little")
        self.N.fromfile(bgen, 4)
        self.N = struct.unpack("<L", self.N)[0]

        if self.N != N:
            raise ValueError('There are {} ids in identifier block != {} ids in headers!'.format(self.N, N))
        if self.identifier_block_len + Lh > offset:
            raise ValueError('Identifier block  + header > Offset !')

        self.ids = OrderedDict()
        for i in range(self.N):
            len_id = ba.bitarray(endian="little")
            len_id.fromfile(bgen, 2)
            len_id = int(len_id.to01()[::-1], 2)
            ids = ba.bitarray(endian="little")
            ids.fromfile(bgen, len_id)
            ids = ids.tostring()
            self.ids[i] = ids


class Bgen_probe(object):
    def __init__(self, bgen, compression, N_ind, layout=None):
        self.compression = compression
        self.layout = layout
        self.N_ind = N_ind
        self.alleles = OrderedDict()
        self.type_dict = {
                    8:np.uint8,
                    16: np.uint16,
                    32: np.uint32,
                    64: np.uint64
                            }

        if self.layout == 1:
            self.N = ba.bitarray(endian="little")
            self.N.fromfile(bgen, 4)
            self.N = struct.unpack("<L", self.N)[0]
            if self.N != self.N_ind:
                raise ValueError('probe #subjects {} != bgen file #{} '.format(self.N, self.N_ind))

        self.Lid = ba.bitarray(endian="little")
        self.Lid.fromfile(bgen, 2)
        self.Lid = struct.unpack("H", self.Lid)[0]
        if self.Lid != 0:
            self.probe_iden = ba.bitarray(endian="little")
            self.probe_iden.fromfile(bgen, self.Lid)
            self.probe_iden = self.probe_iden.tostring()
        else:
            self.probe_iden = 'None'

        self.Len_rsid = ba.bitarray(endian="little")
        self.Len_rsid.fromfile(bgen, 2)
        self.Len_rsid = int(self.Len_rsid.to01()[::-1], 2)

        self.rsid = ba.bitarray(endian="little")
        self.rsid.fromfile(bgen, self.Len_rsid)
        self.rsid = self.rsid.tostring()

        self.len_chr = ba.bitarray(endian="little")
        self.len_chr.fromfile(bgen, 2)
        self.len_chr = int(self.len_chr.to01()[::-1], 2)

        self.CHR = ba.bitarray(endian="little")
        self.CHR.fromfile(bgen, self.len_chr)
        self.CHR = int(self.CHR.tostring())

        self.POS = ba.bitarray(endian="little")
        self.POS.fromfile(bgen, 4)
        self.POS = int(self.POS.to01()[::-1], 2)

        if self.layout == 2:
            self.K = ba.bitarray(endian="little")
            self.K.fromfile(bgen, 2)
            self.K = int(self.K.to01()[::-1], 2)

            for i in range(1, self.K + 1):
                len_a = ba.bitarray(endian="little")
                len_a.fromfile(bgen, 4)
                len_a = int(len_a.to01()[::-1], 2)
                a = ba.bitarray(endian="little")
                a.fromfile(bgen, len_a)
                a = a.tostring()
                self.alleles[i] = [a]

            self.genotypes_length = ba.bitarray(endian="little")
            self.genotypes_length.fromfile(bgen, 4)
            self.genotypes_length = int(self.genotypes_length.to01()[::-1], 2)

            if self.compression != 'raw':
                self.genotypes_length_uncom = ba.bitarray(endian="little")
                self.genotypes_length_uncom.fromfile(bgen, 4)
                self.genotypes_length_uncom = int(self.genotypes_length_uncom.to01()[::-1], 2)

            self.genotypes_bit = ba.bitarray(endian="little")
            need2read = self.genotypes_length if self.compression == 'raw' else self.genotypes_length - 4
            try:
                self.genotypes_bit.fromfile(bgen, need2read)
                self.genotypes_bit = self.genotypes_bit.tobytes()
            except Exception, e:
                print e

        if self.layout == 1:
            self.K = 2
            for i in range(1, self.K + 1):
                len_a = ba.bitarray(endian="little")
                len_a.fromfile(bgen, 4)
                len_a = int(len_a.to01()[::-1], 2)
                a = ba.bitarray(endian="little")
                a.fromfile(bgen, len_a)
                a = a.tostring()
                self.alleles[i] = [a]

            if self.compression == 'raw':
                self.genotypes_length = 6 * self.N
            else:
                self.genotypes_length = ba.bitarray(endian="little")
                self.genotypes_length.fromfile(bgen, 4)
                self.genotypes_length = int(self.genotypes_length.to01()[::-1], 2)

            self.genotypes_bit = ba.bitarray(endian="little")
            try:
                self.genotypes_bit.fromfile(bgen, self.genotypes_length)
                self.genotypes_bit = self.genotypes_bit.tobytes()
            except Exception, e:
                print e

    def info(self):
        print (
        "Iden: {}, RSID: {}, CHR: {}, POS: {}, Alleles: {}".format(self.probe_iden, self.rsid, self.CHR, self.POS,
                                                                   self.alleles))

    def decompression(self):
        if self.genotypes_bit is None:
            raise ValueError('Genotypes data not read yet!')
        if self.compression is None:
            raise ValueError('Compression flag not defined')

        if self.compression == 'zlib':
            return zlib.decompress(self.genotypes_bit)
        elif self.compression == 'zstd':
            return zstd.decompress(self.genotypes_bit)
        elif self.compression == 'raw':
            return self.genotypes_bit


    def get_genotypes(self, ploidy=False, genotypes=False):

        self.genotypes = []
        self.prob = []

        if self.layout == 1:
            genotypes_byte_unzip = self.decompression()
            if len(genotypes_byte_unzip) != self.genotypes_length_uncom:
                raise ValueError('Uncompressed length != length from headers')
            genotypes_unzip_splited = [genotypes_byte_unzip[i:i + 2] for i in range(0, len(genotypes_byte_unzip), 2)]
            for i in genotypes_unzip_splited:
                self.genotypes.append(struct.unpack("H", i)[0] / 32768.)
            self.genotypes = np.array(self.genotypes, dtype=float)
            self.genotypes = np.array(np.split(self.genotypes, 3))
            self.genotypes = (np.array([2, 1, 0]).reshape(3, 1) * self.genotypes).sum(axis=0)
        else:
            genotypes_byte_unzip = self.decompression()
            if len(genotypes_byte_unzip) != self.genotypes_length_uncom:
                raise ValueError('Uncompressed length != length from headers')
            self.N = struct.unpack("<L", genotypes_byte_unzip[:4])[0]
            if self.N != self.N_ind:
                raise ValueError('probe #subjects {} != bgen file #{} '.format(self.N, self.N_ind))
            Ktest = struct.unpack("H", genotypes_byte_unzip[4:6])[0]
            if self.K != Ktest:
                raise ValueError('Alleles #{} !=  alleles #{} in layout '.format(self.K, Ktest))

            self.Pmin = struct.unpack("B", genotypes_byte_unzip[6])[0]
            self.Pmax = struct.unpack("B", genotypes_byte_unzip[7])[0]

            if ploidy:
                self.ploidy = np.fromstring(genotypes_byte_unzip[8:8 + self.N], dtype=np.uint8)
                self.ploidy[self.ploidy>=2**7]=0

                if np.logical_and((self.ploidy != 0), (self.ploidy != 2)).any():
                    raise ValueError('Such ploidy is not supported {}!'.format(self.ploidy))

            self.phased = struct.unpack("B", genotypes_byte_unzip[7 + self.N + 1])[0]

            self.B = ba.bitarray()
            self.B.fromstring(genotypes_byte_unzip[7 + self.N + 1 + 1])
            self.B = int(self.B.to01(), 2)

            if self.phased == 1:
                raise ValueError('Phased data read is not implemented!')
            else:
                prob_unzip = genotypes_byte_unzip[7 + self.N + 1 + 1 + 1:]
                if self.B in self.type_dict:
                    self.prob = np.fromstring(prob_unzip, dtype=self.type_dict[self.B])
                else:
                    if self.B <= 20:
                        decoder = {i: ba.bitarray(format(i, '0{}b'.format(self.B))) for i in range(2 ** self.B)}
                        prob_unzip_splited = ba.bitarray(endian="little")
                        prob_unzip_splited.frombytes(prob_unzip)
                        self.prob = np.array(prob_unzip_splited.decode(decoder), dtype="int")

                self.C = float(2 ** self.B - 1)
                self.prob = self.prob / self.C

                if genotypes:
                    if self.prob.shape[0] / self.N != 2:
                        raise ValueError('Should be two probabilities per subject!')
                    M = self.prob.reshape(self.N, -1).T
                    self.genotypes = (np.array([2, 1]).reshape(2, 1) * M).sum(axis=0)


class Bgen(object):

    def __init__(self, path):
        self.path = path
        self.name = os.path.basename(self.path)
        self.bgen = None
        self.bgen = open(self.path)
        self.size = os.path.getsize(self.path)
        print ('File zise is {} bytes'.format(self.size))
        self.identifier_block =None
        self._indices=False
        self.bgi_txt=None

        try:
            self.offset = ba.bitarray(endian="little")
            self.offset.fromfile(self.bgen, 4)
            self.offset = struct.unpack("<L", self.offset)[0]

            self.header_l = ba.bitarray(endian="little")
            self.header_l.fromfile(self.bgen, 4)
            self.header_l = struct.unpack("<L", self.header_l)[0]

            self.N_probes = ba.bitarray(endian="little")
            self.N_probes.fromfile(self.bgen, 4)
            self.N_probes = struct.unpack("<L", self.N_probes)[0]
            print ('There are {} variants'.format(self.N_probes))

            self.N_ind = ba.bitarray(endian="little")
            self.N_ind.fromfile(self.bgen, 4)
            self.N_ind = struct.unpack("<L", self.N_ind)[0]
            print ('There are {} individuals'.format(self.N_ind))

            self.Magic_number = ba.bitarray(endian="little")
            self.Magic_number.fromfile(self.bgen, 4)
            self.Magic_number = self.Magic_number.tostring()
            if self.Magic_number != 'bgen' and self.Magic_number != '000':
                raise ValueError('Wrong magic number {}'.format(self.Magic_number))

            if (self.header_l - 20) != 0:
                self.Free_data = ba.bitarray(endian="little")
                self.Free_data.fromfile(self.bgen, self.header_l - 20)
            else:
                self.Free_data = None

            self.flags = ba.bitarray(endian="little")
            self.flags.fromfile(self.bgen, 4)
            self.flags = self.flags.to01()

            if int(self.flags[0:2][::-1], 2) == 0:
                self.compression = 'raw'
            elif int(self.flags[0:2][::-1], 2) == 1:
                self.compression = 'zlib'
            elif int(self.flags[0:2][::-1], 2) == 2:
                self.compression = 'zstd'

            layout = int(self.flags[2:6][::-1], 2)
            if layout == 0:
                raise ValueError('layout==0 not supported!')
            elif layout == 1:
                self.layout = 1
            elif layout == 2:
                self.layout = 2
            else:
                raise ValueError('layout > 2 -{}- not supported!'.format(layout))
            print ('Genotype block layout {}'.format(self.layout))

            if self.flags[31] == '0':
                self.identifier_block_len = None
            else:
                self.identifier_block = Bgen_identifier(self.bgen, self.N_ind, self.header_l, self.offset)

            self.start_probes = self.offset + 4
            self.iter_pointer = self.start_probes
            self.bgen.close()
            self.probes_info = OrderedDict()
            self.rsid_info={}

        except Exception, e:
            print e
            self.bgen.close()

    def seek(self, f, pointer):
        if self.size <= pointer:
            raise ValueError('Seek outsize the file!')
        else:
            f.seek(pointer, 0)

    def get_next_probe(self):
        if len(self.probes_info) == self.N_probes:
            return None
        with open(self.path) as f:
            self.seek(f, self.iter_pointer)
            probe = Bgen_probe(f, self.compression, self.N_ind, layout=self.layout)
            self.probes_info[probe.probe_iden ] = [self.iter_pointer, f.tell()]
            if self.rsid_info.get(probe.rsid, 0) == 0:
                self.rsid_info[probe.rsid]=[probe.probe_iden]
            else:
                self.rsid_info[probe.rsid].append(probe.probe_iden)
            self.iter_pointer = f.tell()
        return probe

    def read_probe(self, rsid=None, start=None, iden=None):
        with open(self.path) as f:
            if iden is None:
                start = self.probes_info.get(iden, None)[0]
                if start is None:
                    raise ValueError('IDEN is not known!')
            elif rsid is not None:
                iden=self.rsid_info.get(rsid, None)
                if iden is None:
                    raise ValueError('RSID is not known!')
                elif len(iden)>1:
                    print ('There are several RSID for such IDEN, return all of them!')
                    _=[]
                    for i in iden:
                        self.seek(f, self.probes_info.get(i, None)[0])
                        _.append(    Bgen_probe(f, self.compression, self.N_ind, layout=self.layout)       )
                    return _
                else:
                    self.seek(f, self.probes_info.get(iden, None)[0])
                    probe = Bgen_probe(f, self.compression, self.N_ind, layout=self.layout)
            if start is not None:
                self.seek(f, start)
                probe = Bgen_probe(f, self.compression, self.N_ind, layout=self.layout)

        return probe

    def read_probe_bgi(self,chr=None, pos=None, rsid=None):
        if self.bgi_txt is not None:
            with open(self.path) as f:
                if rsid is not None:
                    if chr is not None and pos is not None:
                        probes = self.bgi_txt.query('chromosome=={} and position=={} and rsid=="{}"'.format(chr, pos, rsid))
                    else:
                        probes=self.bgi_txt.query('rsid=="{}"'.format(rsid))
                elif chr is not None and pos is not None:
                    probes = self.bgi_txt.query('chromosome=={} and position=={}'.format(chr,pos))

                if probes.shape[0]>1:
                    print ('There are several RSID, return all of them!')
                    _ = []
                    for i in probes.file_start_position:
                        self.seek(f, i)
                        _.append(Bgen_probe(f, self.compression, self.N_ind, layout=self.layout))
                    return _
                elif  probes.shape[0]==1:
                    self.seek(f, probes.file_start_position.tolist()[0])
                    probe = Bgen_probe(f, self.compression, self.N_ind, layout=self.layout)
                    return probe
                else:
                    print ('Can not find variant!')
                    return None
        else:
            raise ValueError('bgi file is not defined!')

    def load_bgi_txt(self, path):
        self.bgi_txt=pd.read_csv(path,sep=',')

    def get_indices(self):
        if not self._indices:
            while True:
                try:
                    probe_tmp=self.get_next_probe()
                    if probe_tmp is None:
                        self._indices=True
                        break
                except Exception, e:
                    self._indices = True
                    print e
                    break

    def save_indices(self, path):
        if self._indices:
            dict_tmp={'name':self.name,'probes':self.probes_info, 'rsid':self.rsid_info}
            np.save(os.path.join(path,self.name+ '_ind.npy'),dict_tmp)


    def load_indices(self,path, ovewrite=False):
        dict_tmp=np.load(path).item()
        if self.name!=dict_tmp['name']:
            raise ValueError("Indices file name {} != bgen {} name ".format(dict_tmp['name'], self.name))

        else:
            if len(self.probes_info)!=0 and not ovewrite:
                raise ValueError('Porbes info not empty! If you want to load indices anyway use ovewrite=True parameter.')

            else:
                self.probes_info=dict_tmp['probes']
                self.rsid_info = dict_tmp['rsid']


    def info(self):
        print ("Name:{}; N samples:{}; N probes:{}; Compression:{}; Layout:{}".format(self.name,self.N_ind, self.N_probes, self.compression, self.layout))

