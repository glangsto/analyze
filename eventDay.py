# python file to summarize a directory
# HISTORY
# 25Dec24 GIL Clean up printing of FITS to events
# 25Dec22 GIL add log of fits to event groups
# 25Dec15 GIL fix references to .pdfs
# 25Dec15 GIL initial version

import shutil
import os
from astropy.time import Time
import utcOffset

def htmlFits( html_content, nFit, fits, mjdRef):
    """
    htmlFits() adds successful fits to groups of events to the event log
    where
    html_content  = existing string of html content.
                    output will be appended to this string
    nFit     = number of succesful fits
    fits     = array of fit objects
    mjdRef   = reference MJD for this observing session
    """
    # start showing text without formatting
    if nFit > 0:
        html_content += f"""<pre>\n"""
        if nFit == 1:
            html_content += f"""Successfully fit 1 Event Group.\n"""
        else:
            html_content += f"""Successfully fit {nFit} Event Groups.\n"""
    else:
        html_content += f"""<br>No Fits were made to Event Groups.\n"""
        return html_content

    html_content += "When many flashes are detected we fit the average time of events\n"
    html_content += "FWHM is the time interval of gaussain fit, in hours\n"
    fitStr = "#     Date        Fit Hour    FWYM  HHmmSS.sss         MJD        Event Counts\n"
    html_content += fitStr
    for iFit in range(nFit):
        time = fits[iFit]['time']
        mjd = mjdRef + time

        # Convert MJD to astropy Time object (assumes UTC scale)
        time_obj = Time(mjd, format='mjd', scale='utc')

        # Convert astropy Time object to Python datetime object (in UTC)
        utc_datetime = time_obj.to_datetime()

        # Convert datetime object to a floating-point Unix timestamp (seconds since epoch)
        utc_float_timestamp = utc_datetime.timestamp()
        utcYYmmmDD = utcOffset.strAst( utc_datetime)

        utcstr = str(utc_datetime)

        peak  = fits[iFit]['peak']
        rms   = fits[iFit]['rms']
        stdDev  = fits[iFit]['stdDev']
        fwhm  = fits[iFit]['fwhm']
        utcYYparts = utcYYmmmDD.split('T')
        fitStr = "#FIT %s %7.3f+/-%.3f %6.3f %s     %12.4f %7.1f+/-%4.1f \n" % \
              (utcYYparts[0], time*24., stdDev*24.,\
               fwhm*24., utcYYparts[1], mjd, peak, rms)
        html_content += fitStr
        html_content += "\n"
    html_content += f"""</pre>\n"""  
    return html_content
    # end of htmlFits()

def htmlGroups( htmlFile, calendar, nDir, \
                nGroup, nTen, nHundred, nThousand, \
                nFit, fits, mjdRef):
    """
    htmlGroups() creates a web page for this days events
    where
    htmlFile= open file add the HTML text
    calendar= calendar day of observation (ie 25Dec22)
    nGroup  = total number of groups found
    nTen    = number of groups with between 10 and 99 members   (X)
    nHundred = number of groups with 100 to 999 members         (C)
    nThousand = number of groups with more than 999 members     (M)
    n
    """

    if nGroup <= 0:
        print("<p><h2>No Events matched on %s</h2>" % (calendar),
          file=htmlFile)
    else:
        print("<h2>Summary of events detected by %d telescopes on %s<h1>" %
          (nDir, calendar), file=htmlFile)

        print("<p><table border='2'>",
              file=htmlFile)
        print("<tr><td>Day</td><td>Isolated</td><td>Few</td><td>Groups</td><td>Groups</td><td>Groups of</td></tr>",
              file=htmlFile)
        print("<tr><td>    </td><td>Events</td>  <td>Events</td><td>10 or more</td><td>100 or more</td><td> 1000 or more</td></tr>",
              file=htmlFile)

        print("<tr><td> %s</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td></tr>" % (calendar, nGroup, nTen, nHundred, nThousand),
              file=htmlFile)
        # finish table

    print( "<br>", file=htmlFile)
    
    return
    # end of htmlGroups()

def eventDaySummary(directory_path, calendar, nAll,
                    nFit, fits, mjdRef, output_filename="index.html"):
    """
    Reads directory contents and generates an HTML summary file.
    where directory_path is the location of the directory to summarize
    calendar = calendar day corresponding to events
    nAll = number of events found
    output_filename = summary html file
    nFit = number of successful fits to groups
    fits = fir parameters, a multidimensional list
    mjdRef = reference mjd of fit.
    """
    # set style for plots    
    mystyle = 'style="width: 50%; height: auto;"'
    # Start HTML content
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Event Detected on: {calendar}</title>
    <style>
        body {{ font-family: sans-serif; }}
        ul {{ list-style-type: none; padding-left: 20px; }}
        li.dir {{ font-weight: bold; margin-top: 10px; }}
        li.file {{ margin-left: 20px; }}
    </style>
</head>
<body>
    <h2>On {calendar}: <code>{nAll}</code> Events were detected</h2>
    <p>
    <a href="match-{calendar}.pdf"> <img src="match-{calendar}.svg" {mystyle} alt="Match Histogram"> </a>
    <p>
    <ul>
"""
    # Use os.walk to traverse the directory and its subdirectories
    # The topdown=True argument ensures that we visit a directory before its subdirectories
    for dirpath, dirnames, filenames in os.walk(directory_path):
        # Determine the current level of indentation for HTML
        level = dirpath.replace(directory_path, '').count(os.sep)
        indent = '&nbsp;&nbsp;&nbsp;&nbsp;' * level

        # Add currenthr directory as a list item
        dirparts = dirpath.split("-events-")
        if len(dirparts) == 2:     # if an events directory
            previousDir = os.path.basename(dirpath) or dirpath
            html_content += f'<li class="dir">{indent}📁 <a href="{os.path.basename(dirpath) or dirpath}"> **{os.path.basename(dirpath) or dirpath}**</a> </li>\n'
        else:
            html_content += f'<li class="dir">{indent}📁 **{os.path.basename(dirpath) or dirpath}**</li>\n'
                           
        # Add files in the current directory as list items
        for filename in filenames:
            fileparts = filename.split(".")
            nparts = len(fileparts)
            dirparts = filename.split("/")
            ndir = len(dirparts)
            if nparts > 1:
                if fileparts[1] == 'png':
                    newhref = fileparts[0] + ".pdf"
                    html_content += f'<li class="file">{indent}&nbsp;&nbsp;&nbsp;📄<a href="{newhref}"> <img src="{filename}" {mystyle} alt="{fileparts[0]}"</a></li>\n'
                elif fileparts[1] == 'pdf':
                    continue
                elif fileparts[1] == 'txt':
                    newhref=previousDir + "/" + filename
                    html_content += f'<li class="file">{indent}&nbsp;&nbsp;&nbsp;📄<a href="{newhref}"> {fileparts[0]}</a></li>\n'
                else:
                    continue
    # End HTML content

    more_content = htmlFits( html_content, nFit, fits, mjdRef)

    more_content += '<p><a href="https://www.gb.nrao.edu/~glangsto/flashes/index.html"> Overview of Radio Detection of Cosmic Ray Flashes</a> and access to observations.</p>'

    # Write the content to an HTML file
    with open(output_filename, 'w', encoding='utf-8') as f:
        f.write(more_content)

    f.close()

    outparts = output_filename.split("index.html")
    calendarName = outparts[0] + calendar + ".html"
    print("Copying html file to %s" % (calendarName))
    with open(calendarName, 'w', encoding='utf-8') as ff:
        ff.write(more_content)

    ff.close()

    print(f"Successfully generated HTML summary at {os.path.abspath(output_filename)}")
    print(f"Successfully generated HTML summary at {os.path.abspath(calendarName)}")

    # end of eventDaySummary()

def main():
    """
    Main executable for matching transient events
    """
    inDirName = '/home/karl/gdrive/match/2025/25Jan/25Jan01'
    calendar = '25Jan01'
    mjdRef = 60680.
    nGroup = 1
    outSummaryName = inDirName + "/index.html"
    nFit = 1
    fits = []
    aFit = { 'time': 0.5, 'peak': 1.0, 'stdDev': 0.10, 'fwhm': 1.0, 'rms': 0.1, 'goodFit': True }

    fits.append(aFit)
    
    eventDaySummary('/home/karl/match/2025/25Jan/25Jan01', calendar, nGroup,
                    nFit, fits, mjdRef, outSummaryName)
    
    
if __name__ == "__main__":
    main()

