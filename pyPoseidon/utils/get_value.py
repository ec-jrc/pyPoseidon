"""
Utility function

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon, a software written by George Breyiannis (JRC E.1)
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 

# 
def get_value(mod,kwargs,key,default):
        if key in kwargs.keys():
            return kwargs.get(key)
        elif key in mod.__dict__.keys():
            return mod.__dict__[key]
        else: 
            return default       

