#!/usr/bin/env python

# Convert DES_BCC Galaxy catalogs (Risa Wechler et al.) to a Root tree

import os
import argparse

import numpy as np
from astropy.io import fits
from root_numpy import array2root

parser= argparse.ArgumentParser(description="Convert a DES_BCC Galay Catalog (fits) into a Root tree")
parser.add_argument("--input", action="store", help="input file path")
args = parser.parse_args()

ext = os.path.splitext(args.input)[1]
fn = os.path.splitext(args.input)[0].split("/")[-1]
path = os.path.dirname(args.input)
output = os.path.join(path,fn + ".root")

hdulist = fits.open(args.input)

bcc = hdulist[1].data
length = len(bcc)

# Root N-Tuple content description
root = np.zeros(length, dtype=[('id','i4'),('index','i4'),('ra','f4'),('dec','f4'),('z','f4'),('gamma1','f4'),('gamma2','f4'),('kappa','f4'),
    ('size','f4'),('eps1','f4'),('eps2','f4'),('mag','f4'),('teps1','f4'),('teps2','f4'),('tra','f4'),('tdec','f4'),('mu','f4'),('tsize','f4')])


root['id'] = bcc['ID']
print "ID done..."
root['index'] = bcc['INDEX']
print "INDEX done..."
root['ra'] = bcc['RA']
print "RA done..."
root['dec'] = bcc['DEC']
print "DEC done"
root['z'] = bcc['Z']
print "Z done"
root['gamma1'] = bcc['GAMMA1']
print "GAMMA1 done..."
root['gamma2'] = bcc['GAMMA2']
print "GAMMA2 done..."
root['kappa'] = bcc['KAPPA']
print "KAPPA done..."
root['size'] = bcc['SIZE']
print "SIZE done..."
root['eps1'] = bcc["EPSILON"][0:,0]
print "EPSILON 1 done..."
root['eps2'] = bcc["EPSILON"][0:,1]
print "EPSILON 2 done..."
root["mag"] = bcc["TMAG"][0:,2]
print "TMAG done..."
root["teps1"] = bcc["TE"][0:,0]
print "TEPS1 done..."
root["teps2"] =bcc["TE"][0:,1]
print "TEPS2 done..."
root["tra"] = bcc["TRA"]
print "TRA done..."
root["tdec"] = bcc["TDEC"]
print "TDEC done..."
root["tsize"] = bcc["TSIZE"]
print "TSIZE done..."
root["mu"] = bcc["MU"]

print "All Done !"

array2root(root,output,'bcc')




