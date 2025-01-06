#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

from bronx.syntax.decorators import nicedeco


@nicedeco
def init_before(mtd):  # TODO: remove when integrated in bronx (PR submitted)
    """
    Decorator for methods: call method self._actual_init()
    before actually calling method if not self.initialized.
    """
    def initialized(self, *args, **kwargs):
        if not self._initialized:
            self._actual_init()
        return mtd(self, *args, **kwargs)
    return initialized


