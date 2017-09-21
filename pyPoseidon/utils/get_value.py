# 
def get_value(mod,kwargs,key,default):
        try:
            return kwargs.get(key, mod.__dict__[key])
        except KeyError: 
            return default       
