# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 10:36:35 2020

@author: sp3660
"""

def PathStructure():
    import os
    
    global ProjectCodePath
    ProjectCodePath=os.getcwd()
    global ProjectName
    ProjectName=os.path.basename(ProjectCodePath)
    global ProjectDataPath
    ProjectDataPath='G:\\Dropbox\\CodeData\\'+ ProjectName
    global ProjectRAWPath
    ProjectRAWPath='F:\\CodeRawData\\'+ ProjectName
    global ProjectTempRAWPath
    ProjectTempRAWPath='G:\\CodeTempRawData\\'+ ProjectName
    return ProjectCodePath, ProjectName, ProjectDataPath, ProjectRAWPath, ProjectTempRAWPath