#!/usr/bin/env python3
import os

def createdir(relative_path):
    if not os.path.exists(relative_path):
        os.mkdir(relative_path)

