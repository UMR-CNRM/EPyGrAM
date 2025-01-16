#!/usr/bin/env python
# -*- coding: utf-8 -*
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
This module aims at reading tiff with private IFDs.
It uses PIL for image reading.

This module uses code from pylibtiff
(https://pypi.python.org/pypi/libtiff, https://code.google.com/p/pylibtiff or https://github.com/hmeine/pylibtiff)
"""

import os
import numpy
import mmap
import PIL.Image
import io


class PyexttiffError(Exception):

    """Error handling for pyexttiff."""
    pass


class TiffFile(object):

    """
    This class allows the access to the entire tiff file (tags and images).
    """

    #see https://github.com/numpy/numpy/issues/2407 for the str workaround
    #We will be able to suppress it when python2 with "from __future__ import unicode_literals" will not be used
    _rational = numpy.dtype([(str('numer'), numpy.uint32), (str('denom'), numpy.uint32)])
    _srational = numpy.dtype([(str('numer'), numpy.int32), (str('denom'), numpy.int32)])

    _type2name = {1:'BYTE', 2:'ASCII', 3:'SHORT', 4:'LONG', 5:'RATIONAL',  # two longs, lsm uses it for float64
                  6:'SBYTE', 7:'UNDEFINED', 8:'SSHORT', 9:'SLONG', 10:'SRATIONAL',
                  11:'FLOAT', 12:'DOUBLE',
                  }
    _name2type = dict((v, k) for k, v in _type2name.items())
    _name2type['SHORT|LONG'] = _name2type['LONG']
    _name2type['LONG|SHORT'] = _name2type['LONG']
    _type2bytes = {1:1, 2:1, 3:2, 4:4, 5:8, 6:1, 7:1, 8:2, 9:4, 10:8, 11:4, 12:8}
    _type2dtype = {1:numpy.uint8, 2:numpy.uint8, 3:numpy.uint16, 4:numpy.uint32, 5:_rational,
                   6:numpy.int8, 7:numpy.uint8, 8:numpy.int16, 9:numpy.int32, 10:_srational,
                   11:numpy.float32, 12:numpy.float64}

    class _LittleEndianNumpyDTypes(object):
        uint8 = numpy.dtype('<u1')
        uint16 = numpy.dtype('<u2')
        uint32 = numpy.dtype('<u4')
        uint64 = numpy.dtype('<u8')
        int8 = numpy.dtype('<i1')
        int16 = numpy.dtype('<i2')
        int32 = numpy.dtype('<i4')
        int64 = numpy.dtype('<i8')
        float32 = numpy.dtype('<f4')
        float64 = numpy.dtype('<f8')
        complex64 = numpy.dtype('<c8')
        complex128 = numpy.dtype('<c16')

        @property
        def type2dt(self):
            return dict((k, numpy.dtype(v).newbyteorder('<')) for k, v in TiffFile._type2dtype.items())

    class _BigEndianNumpyDTypes(object):
        uint8 = numpy.dtype('>u1')
        uint16 = numpy.dtype('>u2')
        uint32 = numpy.dtype('>u4')
        uint64 = numpy.dtype('>u8')
        int8 = numpy.dtype('>i1')
        int16 = numpy.dtype('>i2')
        int32 = numpy.dtype('>i4')
        int64 = numpy.dtype('>i8')
        float32 = numpy.dtype('>f4')
        float64 = numpy.dtype('>f8')
        complex64 = numpy.dtype('>c8')
        complex128 = numpy.dtype('>c16')

        @property
        def type2dt(self):
            return dict((k, numpy.dtype(v).newbyteorder('>')) for k, v in TiffFile._type2dtype.items())

    def __init__(self, filename, subIFDpaths=[], method=1):
        """
        Opens a tiff file, reads header and IFDs.
        *filename* is the filename containing the tiff
        *subIFDpaths* is the list of tag path whose values are offset to private IFDs
            a tag path is a tuple representing the path to a given tag which must represent an IFD
            (34665) means that tag 34665 of any given public IFD is an offset to a private IFD
            (32001, 521) means that tag 32001 of any given public IFD is an offset to a private IFD
                           and that tag 521 of any private tag referenced by a 32001 public tag is also an offset to a private IFD
        *method* is the method used to read data:
            1: f=open(..., 'rb') ; numpy.frombuffer(f.read(), dtype=numpy.ubyte)
            2: f=open(..., 'rb') ; numpy.ndarray(buffer=mmap(f), dtype=numpy.ubyte)
            3: same as 2 but with modifications allowed - DANGEROUS
        """

        self._filename = filename
        self._subIFDpaths = subIFDpaths
        self._fileHandle = None
        if not os.path.exists(filename):
            raise IOError(filename + " must exists.")

        # Reading file
        if method == 1:
            with open(filename, 'rb') as f:
                self._data = numpy.frombuffer(f.read(), dtype=numpy.ubyte)
        elif method == 2:
            self._fileHandle = open(filename, 'rb')
            mm = mmap.mmap(self._fileHandle.fileno(), 0, prot=mmap.PROT_READ)
            self._data = numpy.ndarray(shape=(mm.size(),), buffer=mm, dtype=numpy.ubyte)
        elif method == 3:
            self._fileHandle = open(filename, 'r+b')
            mm = mmap.mmap(self._fileHandle.fileno(), 0)
            self._data = numpy.ndarray(shape=(mm.size(),), buffer=mm, dtype=numpy.ubyte)
        else:
            raise PyexttiffError("This method is unknown.")

        # Decoding header - byte order
        byteorder = self._data[0:2].view(dtype=numpy.uint16)[0]
        if byteorder == 0x4949:
            self.endian = 'little'
            self.dtypes = TiffFile._LittleEndianNumpyDTypes()
        elif byteorder == 0x4d4d:
            self.endian = 'big'
            self.dtypes = TiffFile._BigEndianNumpyDTypes()
        else:
            raise IOError('unrecognized byteorder: %s' % (hex(byteorder)))

        # Decoding header - magic number
        magic = self._get_uint16(2)
        if magic != 42:
            raise IOError('wrong magic number for TIFF file: %s' % (magic))

        # Decoding header - first IFD offset
        IFD0offset = self._get_uint32(4)

        # IFD reading
        self.IFDs = []
        offset = IFD0offset
        num = 0
        while offset:
            ifd, offset = self._readIFD(offset, (), subIFDpaths, num)
            if ifd.has_image():
                num += 1
            else:
                raise PyexttiffError("Not sure about IFD that does not contain image.")
            self.IFDs.append(ifd)

    def _readIFD(self, offset, path, subIFDpaths, num):
        """
        Reads recursively IFDs
        """
        ifd = IFD()
        n = self._get_uint16(offset)
        for i in range(n):
            entryOffset = offset + 2 + i * 12
            entrytag = self._get_uint16(entryOffset)
            entrytype = self._get_uint16(entryOffset + 2)
            entrycount = self._get_uint32(entryOffset + 4)
            entrybytes = TiffFile._type2bytes.get(entrytype, 0)
            if entrycount == 1 and 1 <= entrybytes <= 4:
                entryvalue = self._get_values(entryOffset + 8, entrytype, entrycount)
            else:
                valueOffset = self._get_int32(entryOffset + 8)
                entryvalue = self._get_values(valueOffset, entrytype, entrycount)
            entrypath = tuple(list(path) + [entrytag])
            if entrypath in [mypath[:len(entrypath)] for mypath in subIFDpaths]:
                #In the first version of this tool, entryvalue was an array at this stage
                #This is corrected to suppress the numpy warning but code
                #(by replacing entryvalue by entryvalue[0] in following statement)
                #lacks comment and I'm now unable to understand this part.
                #If an error is raised by this line, it would be necessary
                #to investigate more...
                assert entryvalue.shape == (1, ), "Not as expected..."
                
                subifd, _ = self._readIFD(entryvalue[0], entrypath, subIFDpaths, None)
                ifd.append(IFDEntry(entrytag, entrytype, subifd))
            else:
                ifd.append(IFDEntry(entrytag, entrytype, entryvalue))
        if path == () and ifd.has_tag(273):
            # Raw data
            nbRows = ifd.get_value(257)
            offsetValues = ifd.get_value(273)
            nbRowsPerStrip = ifd.get_value(278)
            nbBytesPerStrip = ifd.get_value(279)
            if not isinstance(offsetValues, numpy.ndarray):
                offsetValues = numpy.array([offsetValues])
                nbBytesPerStrip = numpy.array([nbBytesPerStrip])
            if nbRows / nbRowsPerStrip + (1 if nbRows % nbRowsPerStrip != 0 else 0) != len(offsetValues):
                raise PyexttiffError("Total number of rows, strip numbers and number of rows per strips are not consistent.")
            data = []
            for i in range(len(offsetValues)):
                data.append(self._get_values(offsetValues[i], 1, nbBytesPerStrip[i]))
            ifd.get_entry(273).set_value(data)

            # Image data
            im = self.get_PILImage()
            im.seek(num)
            data = numpy.array(im)
            if data.shape == ():
                data = numpy.array(im)  # Sometimes must be called twice to return values...
            ifd._image = data
        nextIFD = self._get_uint32(offset + 2 + n * 12)
        return (ifd, nextIFD)

    def get_data(self):
        """
        Returns the ndarray containing the data.
        """
        return self._data

    def get_PILImage(self):
        """
        Returns the PIL image object of the file.
        """
        try:
            meth = numpy.getbuffer
        except AttributeError:
            meth = memoryview
        return PIL.Image.open(io.BytesIO(meth(self.get_data())))

    def _get_uint16(self, offset):
        return self.get_data()[offset:offset + 2].view(dtype=self.dtypes.uint16)[0]
    def _get_uint32(self, offset):
        return self.get_data()[offset:offset + 4].view(dtype=self.dtypes.uint32)[0]
    def _get_int32(self, offset):
        return self.get_data()[offset:offset + 4].view(dtype=self.dtypes.int32)[0]
    def _get_values(self, offset, typ, count):
        if isinstance(typ, numpy.dtype):
            dtype = typ
            size = typ.itemsize
        elif isinstance(typ, type) and issubclass(typ, numpy.generic):
            dtype = typ
            size = typ().itemsize
        else:
            if isinstance(typ, str):
                ntyp = typ
                typ = TiffFile._name2type.get(typ)
            else:
                ntyp = str(typ)
            dtype = self.dtypes.type2dt.get(typ)
            size = TiffFile._type2bytes.get(typ)
            if dtype is None or size is None:
                raise PyexttiffError('_get_values: incomplete info for type=%r [%r]: dtype=%s, bytes=%s\n' % (typ, ntyp, dtype, size))
        result = self.get_data()[offset:offset + size * count].view(dtype=dtype)
        return result

    def close(self):
        """
        Closes the file.
        """
        if self._fileHandle is not None:
            self._fileHandle.close()

    def __del__(self):
        """
        __del__ method
        """
        self.close()


class IFD(list):
    """This class represent an IFD."""

    def __init__(self):
        """Initialisation method of IFD class."""
        self._image = None

    def has_tag(self, tag):
        """
        Returns True if an entry fits the tag given
        :param tag: tag to look for, as an integer or a name
        :return: True if tag exists
        """
        if isinstance(tag, int):
            return tag in self.get_tagValues()
        else:
            return tag in self.get_tagNames()

    def get_entry(self, tag):
        """
        Returns the entry for a tag
        :param tag: tag to look for, as an integer or a name
        :return: the entry associated to the tag
        """
        if not self.has_tag(tag):
            raise PyexttiffError("This tag doesn't exist in this IFD.")
        for entry in self:
            if (entry.get_tagValue() if isinstance(tag, int) else entry.get_tagName()) == tag:
                result = entry
                break
        return result

    def get_value(self, tag, human=True):
        """
        Returns the value for a tag
        :param tag: tag to look for, as an integer or a name
        :param human: if True, value is modified:
                        - value[0] is returned instead of value if array contains only one element
                        - conversion in string is achieved for arrays representing strings
        :return: the value associated to the tag
        """
        if not self.has_tag(tag):
            raise PyexttiffError("This tag doesn't exist in this IFD.")
        return self.get_entry(tag).get_value(human)

    def has_image(self):
        """Returns True if one tag is an image"""
        return self._image is not None

    def get_image(self):
        """Returns the image"""
        if not self.has_image():
            raise PyexttiffError("This IFD doesn't contain an image.")
        return self._image

    def get_tagValues(self):
        """Returns the list of the tags as decimal values"""
        return [entry.get_tagValue() for entry in self]

    def get_tagNames(self):
        """Returns the list of the tag names"""
        return [entry.get_tagName() for entry in self]

    def as_dict(self, keys='value'):
        """
        Returns a dictionary containing all entries.
        :param keys: keys to use for the dictionary, among ('value', 'name')
        :return: the dictionary
        """
        assert keys in ('value', 'name'), "keys must be in ('value', 'name')"
        return {entry.get_tagValue() if keys == 'value' else entry.get_tagName: entry.get_value() for entry in self}

class IFDEntry(object):
    """This class represent an IFD entry"""

    # <TagName> <Hex> <Type> <Number of values>
    _tag_info = '''
# standard tags:
NewSubfileType FE LONG 1
SubfileType FF SHORT 1
ImageWidth 100 SHORT|LONG 1
ImageLength 101 SHORT|LONG 1
BitsPerSample 102 SHORT SamplesPerPixel
Compression 103 SHORT 1
  Uncompressed 1
  CCITT1D 2
  Group3Fax 3
  Group4Fax 4
  LZW 5
  JPEG 6
  PackBits 32773
PhotometricInterpretation 106 SHORT 1
  WhiteIsZero 0
  BlackIsZero 1
  RGB 2
  RGBPalette 3
  TransparencyMask 4
  CMYK 5
  YCbCr 6
  CIELab 8
Threshholding 107 SHORT 1
CellWidth 108 SHORT 1
CellLength 109 SHORT 1
FillOrder 10A SHORT 1
DocumentName 10D ASCII
ImageDescription 10E ASCII
Make 10F ASCII
Model 110 ASCII
StripOffsets 111 SHORT|LONG StripsPerImage
Orientation 112 SHORT 1
  TopLeft 1
  TopRight 2
  BottomRight 3
  BottomLeft 4
  LeftTop 5
  RightTop 6
  RightBottom 7
  LeftBottom 8
SamplesPerPixel 115 SHORT 1
RowsPerStrip 116 SHORT|LONG 1
StripByteCounts 117 LONG|SHORT StripsPerImage
MinSampleValue 118 SHORT SamplesPerPixel
MaxSampleValue 119 SHORT SamplesPerPixel
XResolution 11A RATIONAL 1
YResolution 11B RATIONAL 1
PlanarConfiguration 11C SHORT 1
  Chunky 1
  Planar 2
PageName 11D ASCII
XPosition 11E DOUBLE
YPosition 11F DOUBLE
FreeOffsets 120 LONG
FreeByteCounts 121 LONG
GrayResponseUnit 122 SHORT 1
GrayResponseCurve 123 SHORT 2**BitsPerSample
T4Options 124 LONG 1
T6Options 125 LONG 1
ResolutionUnit 128 SHORT 1
PageNumber 129 SHORT 2
TransferFunction 12D SHORT (1|SamplesPerPixel)*2**BitsPerSample
Software 131 ASCII
DateTime 132 ASCII 20
Artist 13B ASCII
HostComputer 13C ASCII
Predictor 13D SHORT 1
WhitePoint 13E RATIONAL 2
PrimaryChromaticities 13F RATIONAL 6
ColorMap 140 SHORT 3*(2**BitsPerSample)
HalftoneHints 141 SHORT 2
TileWidth 142 SHORT|LONG 1
TileLength 143 SHORT|LONG 1
TileOffsets 144 LONG TilesPerImage
TileByteCounts 145 SHORT|LONG TilesPerImage
InkSet 14C SHORT 1
InkNames 14D ASCII <total number of chars in all ink name strings, including zeros>
NumberOfInks 14E SHORT 1
DotRange 150 BYTE|SHORT 2|2*NumberOfInks
TargetPrinter 151 ASCII any
ExtraSamples 152 BYTE <number of extra components per pixel>
SampleFormat 153 SHORT SamplesPerPixel
SMinSampleValue 154 Any SamplesPerPixel
SMaxSampleValue 155 Any SamplesPerPixel
TransferRange 156 SHORT 6
JPEGProc 200 SHORT 1
JPEGInterchangeFormat 201 LONG 1
JPEGInterchangeFormatLength 202 LONG 1
JPEGRestartInterval 203 SHORT 1
JPEGLosslessPredictos 205 SHORT SamplesPerPixel
JPEGPointTransforms 206 SHORT SamplesPerPixel
JPEGQTables 207 LONG SamplesPerPixel
JPEGDCTables 208 LONG SamplesPerPixel
JPEGACTables 209 LONG SamplesPerPixel
YCbCrCoefficients 211 RATIONAL 3
YCbCrSubSampling 212 SHORT 2
YCbCrPositioning 213 SHORT 1
ReferenceBlackWhite 214 LONG 2*SamplesPerPixel
Copyright 8298 ASCII Any

# non-standard tags:
CZ_LSMInfo 866C CZ_LSM

# EXIF tags, see http://www.awaresystems.be/imaging/tiff/tifftags/privateifd/exif.html
EXIF_IFDOffset 8769 SHORT 1
EXIF_ExposureTime 829a RATIONAL 1
EXIF_FNumber 829d RATIONAL 1
EXIF_ExposureProgram 8822 SHORT 1
EXIF_SpectralSensitivity 8824 ASCII
EXIF_ISOSpeedRatings 8827 SHORT 1
EXIF_OECF 8828 UNDEFINED
EXIF_ExifVersion 9000 UNDEFINED 4
EXIF_DateTimeOriginal 9003 ASCII
EXIF_DateTimeDigitized 9004 ASCII
EXIF_ComponentsConfiguration 9101 UNDEFINED 4
EXIF_CompressedBitsPerPixel 9102 RATIONAL 1
EXIF_ShutterSpeedValue 9201 SRATIONAL 1
EXIF_ApertureValue 9202 RATIONAL 1
EXIF_BrightnessValue 9203 SRATIONAL 1
EXIF_ExposureBiasValue 9204 SRATIONAL 1
EXIF_MaxApertureValue 9205 RATIONAL 1
EXIF_SubjectDistance 9206 RATIONAL 1
EXIF_MeteringMode 9207 SHORT 1
EXIF_LightSource 9208 SHORT 1
EXIF_Flash 9209 SHORT 1
EXIF_FocalLength 920a RATIONAL 1
EXIF_SubjectArea 9214 SHORT 2|3|4
EXIF_MakerNote 927c UNDEFINED
EXIF_UserComment 9286 UNDEFINED
EXIF_SubsecTime 9290 ASCII
EXIF_SubsecTimeOriginal 9291 ASCII
EXIF_SubsecTimeDigitized 9292 ASCII
EXIF_FlashpixVersion a000 UNDEFINED 4
EXIF_ColorSpace a001 SHORT 1
EXIF_PixelXDimension a002 SHORT!LONG 1
EXIF_PixelYDimension a003 SHORT!LONG 1
EXIF_RelatedSoundFile a004 ASCII 13
EXIF_FlashEnergy a20b RATIONAL 1
EXIF_SpatialFrequencyResponse a20c UNDEFINED
EXIF_FocalPlaneXResolution a20e RATIONAL 1
EXIF_FocalPlaneYResolution a20f RATIONAL 1
EXIF_FocalPlaneResolutionUnit a210 SHORT 1
EXIF_SubjectLocation a214 SHORT 2
EXIF_ExposureIndex a215 RATIONAL 1
EXIF_SensingMethod a217 SHORT 1
EXIF_FileSource a300 UNDEFINED 1
EXIF_SceneType a301 UNDEFINED 1
EXIF_CFAPattern a302 UNDEFINED
EXIF_CustomRendered a401 SHORT 1
EXIF_ExposureMode a402 SHORT 1
EXIF_WhiteBalance a403 SHORT 1
EXIF_DigitalZoomRatio a404 RATIONAL 1
EXIF_FocalLengthIn35mmFilm a405 SHORT 1
EXIF_SceneCaptureType a406 SHORT 1
EXIF_GainControl a407 SHORT 1
EXIF_Contrast a408 SHORT 1
EXIF_Saturation a409 SHORT 1
EXIF_Sharpness a40a SHORT 1
EXIF_DeviceSettingDescription a40b UNDEFINED
EXIF_SubjectDistanceRange a40c SHORT 1
EXIF_ImageUniqueID a420 ASCII 33

'''
    _tag_value2name = {}
    _tag_name2value = {}
    _tag_value2type = {}
    for line in _tag_info.split('\n'):
        if not line or line.startswith('#'):
            continue
        if line[0] == ' ':
            pass
        else:
            n, h, t = line.split()[:3]
            h = eval('0x' + h)
            _tag_value2name[h] = n
            _tag_value2type[h] = t
            _tag_name2value[n] = h

    def __init__(self, tag, entrytype, value=None):
        """
        *tag* is the tag number of the entry
        *entrytype* is the type of the entry
        *value* is the value associated to the tag
        """
        self._tag = tag
        self._type = entrytype
        self._value = value
        self._name = self._tag_value2name.get(tag, 'TAG%s' % (hex(tag),))

    def is_image(self):
        """Returns True if content in an image"""
        return self.get_tagValue() == 273

    def is_IFD(self):
        """Returns True if content is an IFD."""
        return isinstance(self.get_value(), IFD)

    def get_tagValue(self):
        """Returns the tag"""
        return self._tag

    def get_value(self, human=True):
        """
        Returns the value
        if human=True, value is modified:
            - value[0] is returned instead of value if array contains only one element
            - conversion in string is achieved for arrays representing strings
        """
        value = self._value
        if human:
            if len(value) == 1:
                value = value[0]
            if self.get_type() == 2: #ASCII
                value = (b''.join(value.view('|S%s' % (value.nbytes // value.size)))).decode('UTF8')
        return value

    def get_tagName(self):
        """Returns the tag name"""
        return self._name

    def get_type(self):
        """Returns the type of entry."""
        return self._type

    def set_value(self, value):
        """Sets the value of the entry."""
        self._value = value
