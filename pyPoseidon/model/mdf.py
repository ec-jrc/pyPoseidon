"""read/write Delft3D-FLOW *.mdf input files to/from dictionary"""

__version__ = "$Revision: 7870 $"

#  Copyright notice
#   --------------------------------------------------------------------
#   Copyright (C) 2012 Deltares
#       Gerben J. de Boer
#
#       gerben.deboer@deltares.nl
#
#   This library is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this library.  If not, see <http://www.gnu.org/licenses/>.
#   --------------------------------------------------------------------
#
# This tool is part of <a href="http://www.OpenEarth.eu">OpenEarthTools</a>.
# OpenEarthTools is an online collaboration to share and manage data and
# programming tools in an open source, version controlled environment.
# Sign up to recieve regular updates of this function, and to contribute
# your own tools.

# $Id: mdf.py 7870 2012-12-31 13:33:52Z boer_g $
# $Date: 2012-12-31 14:33:52 +0100 (Mon, 31 Dec 2012) $
# $Author: boer_g $
# $Revision: 7870 $
# $HeadURL: https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/OpenEarthTools/openearthtools/io/delft3d/mdf.py $
# $Keywords: $

def _RHS2val_(line):
   """parse 1(!) RHS line value from *.mdf file to a str, '' or float"""

   if ('#' in line):
      dummy, value, dummy = line.split('#',2)
   elif ('[' in line):
      value = '';
   else:
      value = line.strip();
      value = map(float, value.split())
   return value

def read(mdf_file):
   """Read Delft3D-Flow *.mdf file into dictionary.
      >> inp, inp_order = mdf.read('a.mdf')
   Note that all comment lines and keyword formatting are not preserved.
   And mind that inp is a dictionary where the file's keyword order
   is not preserved. The order from the file
   is therefore also be returned in variable inp_order.
   This can be ignored in 2 standard python ways:
      >> inp,_ = mdf.read('a.mdf')
      >> inp   = mdf.read('a.mdf')[0]"""
   
   keywords = {} # empty dictionary
   keyword_order = []
   
   f = open(mdf_file,'r')
   
   for line in f.readlines():
   
      if ('=' in line):
         keyword, value = line.split('=',1)
         keyword = keyword.strip()
         value   = value.strip()
         keywords[keyword] = []
         new     = True
         keyword_order.append(keyword)

      elif not(keyword=='Commnt'):
         value = line.lstrip()
         new   = False
         
      if not(keyword=='Commnt'):
         value = _RHS2val_(value)
         if new:
            keywords[keyword] = value
         else:
            if  type(value) is str:
               keywords[keyword] = keywords[keyword] + value
            else:
               keywords[keyword].append(value[0])
               
   f.close()

   return keywords, keyword_order  

def _val2RHS_(f,keyword,value):
   """parse a list of str, '' or floats to multiple (!, if needed) RHS *.mdf lines"""

   # values that need to be written column-wise rather than row wise (although short row vectors are allowed)
   columnwise_list = ['Thick','Rettis','Rettib','u0','v0','s0','t0','c01','c02','c03','c04','c01'];

   if  type(value) is str:
      if (keyword=='Runtxt'):
         w = 30
         f.write(keyword.ljust(7) + '= #' + str(value[:w]) + '#\n')
         for i in range(w,len(value),w):
            f.write('         #' + value[i:i+w] + '#\n')
      else:
         f.write(keyword.ljust(7) + '= #' + value + '#\n')
   else:
      if keyword in columnwise_list:
         f.write(keyword.ljust(7) + '= ' + str(value[0]) + '\n')
         for val in value[1:]:
            f.write('         ' + str(val) + '\n')
      else:
         f.write(keyword.ljust(7) + '= ' + ' '.join(("%g" % x) for x in value) + '\n')

def write(keywords,mdf_file,**kwargs):
   """Write dictionary to Delft3D-FLOW *.mdf file.
   The keywords are written in ramdom order,
      >> mdf.write(inp, 'b.mdf')
   An keyword 'selection' can be passed containing the 
   desired order (for instance resurned by mdf.read())
      >> inp, order = mdf.read('a.mdf')
      >> mdf.write(inp, 'b.mdf',selection=order)
   This can also be used to write only a subset of keywords to file.
      >> mdf.write(inp,'c.mdf',selection=['Runtxt','MNKmax'])
      
   To ignore a keyword use keyword 'exclude', e.g. to enforce cold start:
      >> mdf.write(inp,'c.mdf',exclude=['Restid']) 
      
   Example: modify time step Dt of collection Delft3D-FLOW simulations:
      >> import mdf, glob, os
      >> mdf_list = glob.glob('*.mdf');
      >> for mdf_file in mdf_list:
      >>    inp, ord = mdf.read(mdf_file)
      >>    inp['Dt'] = [1]
      >>    mdf_base, ext = os.path.splitext(mdf_file)
      >>    mdf.write(inp, mdf_base + '_dt=1' + ext, selection=ord)
      """

   OPT = {}
   OPT['selection'] = []
   OPT['exclude']   = []
      
   for k,v in kwargs.iteritems():
       OPT[k] = v
   
   f = open(mdf_file,'w')
   
   if (OPT['selection'] is None):
      for keyword in keywords:
         if not(keyword in OPT['exclude']):
            value = keywords[keyword]
            _val2RHS_(f,keyword,value)
   else:   
      for keyword in OPT['selection']:
         if not(keyword in OPT['exclude']):
            value = keywords[keyword]
            _val2RHS_(f,keyword,value)
      
   f.close()
