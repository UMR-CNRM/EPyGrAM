#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
Contains defaults for GRIB encoding.
"""

# GRIB2 ------------------------------------------------------------------------
#: GRIB2 key/value defaults, ordered by section
GRIB2_keyvalue = {1:{'tablesVersion':15,
                     'productionStatusOfProcessedData':2,
                     'typeOfProcessedData':2},
                  2:{},
                  3:{'shapeOfTheEarth':6,
                     'iScansNegatively':0,
                     'jScansPositively':0,
                     'jPointsAreConsecutive':0},
                  4:{'typeOfGeneratingProcess':2,
                     'generatingProcessIdentifier':255,
                     # validity
                     'hoursAfterDataCutoff':None,
                     'minutesAfterDataCutoff':None,
                     'indicatorOfUnitOfTimeRange':13,
                     'indicatorOfUnitForTimeRange':13,
                     'indicatorOfUnitForTimeIncrement':13,
                     # second fixed surface
                     'typeOfSecondFixedSurface':255,
                     'scaleFactorOfFirstFixedSurface':0,
                     'scaleFactorOfSecondFixedSurface':0,
                     # simulated satellite imagery
                     'NB':1,
                     'satelliteSeries':0,
                     'satelliteNumber':0,
                     'instrumentType':0,
                     'scaleFactorOfCentralWaveNumber':0},
                  5:{'packingType':'grid_second_order',
                     'bitsPerValue':12}
                  }
#: metadata to embark in field.misc_metadata
GRIB2_metadata_to_embark = ['typeOfGeneratingProcess',
                            'productionStatusOfProcessedData',
                            'typeOfProcessedData',
                            'generatingProcessIdentifier',
                            'typeOfTimeIncrement',
                            'indicatorOfUnitOfTimeRange',
                            'indicatorOfUnitForTimeRange',
                            'indicatorOfUnitForTimeIncrement',
                            'hoursAfterDataCutoff',
                            'minutesAfterDataCutoff',
                            'timeIncrement',
                            'bitsPerValue',
                            'dataRepresentationTemplateNumber',
                            'uvRelativeToGrid',
                            ]

# GRIB1 ------------------------------------------------------------------------
#: GRIB1 sample (from epygram)
GRIB1_sample = 'GRIB1_grid_second_order'
#: packing
GRIB1_packing = {'packingType':'grid_second_order',
                 'bitsPerValue':16}
#: ordering
GRIB1_ordering = {'iScansNegatively':0,
                  'jScansPositively':0,
                  'jPointsAreConsecutive':0}
#: GRIB1 key/value defaults, ordered by section
GRIB1_keyvalue = {'generatingProcessIdentifier':GRIB2_keyvalue[4]['generatingProcessIdentifier'],
                  }
