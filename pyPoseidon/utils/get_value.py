# 
def get_value(mod,kwargs,key,default):
        if key in kwargs.keys():
            return kwargs.get(key)
        elif key in mod.__dict__.keys():
            return mod.__dict__[key]
        else: 
            return default       

