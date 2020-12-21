#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 09:44:58 2020

@author: janbinkowski
"""

from typing import Callable
import multiprocessing as mp

class Process:
    def __init__(self, object_to_process: list, number_of_process: int) -> None:
        self.object_to_process = object_to_process
        self.number_of_process = number_of_process
        self.processed = None
                
    def run_processes(self, func: Callable) -> None:
        with mp.Pool(self.number_of_process) as pool:
            result = pool.map(func, self.object_to_process)
            self.processed = result