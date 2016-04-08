from onvif import ONVIFCamera

def run():

    mycam = ONVIFCamera('192.168.1.10', 8899, 'admin', '', '/etc/onvif/wsdl/')

    # Get Hostname
    resp = mycam.devicemgmt.GetHostname()
    print 'My camera`s hostname: ' + str(resp)

    # Get system date and time
    dt = mycam.devicemgmt.GetSystemDateAndTime()
    tz = dt.TimeZone
    year = dt.UTCDateTime.Date.Year
    hour = dt.UTCDateTime.Time.Hour

if __name__ == "__main__":
    run()