import datetime

def myconverter(o):
    if isinstance(o, datetime.datetime):
        return o.__str__()

