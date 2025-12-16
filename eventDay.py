# python file to summarize a directory
# HISTORY
# 25Dec15 GIL initial version

import os

def eventDaySummary(directory_path, calendar, output_filename="index.html"):
    """
    Reads directory contents and generates an HTML summary file.
    where directory_path is the location of the directory to summarize
    calendar = calendar day corresponding to events
    output_filename = summary html file
    """
    
    # Start HTML content
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Directory Listing: {directory_path}</title>
    <style>
        body {{ font-family: sans-serif; }}
        ul {{ list-style-type: none; padding-left: 20px; }}
        li.dir {{ font-weight: bold; margin-top: 10px; }}
        li.file {{ margin-left: 20px; }}
    </style>
</head>
<body>
    <h1>Contents of <code>{directory_path}</code></h1>
    <ul>
"""
    # Use os.walk to traverse the directory and its subdirectories
    # The topdown=True argument ensures that we visit a directory before its subdirectories
    for dirpath, dirnames, filenames in os.walk(directory_path):
        # Determine the current level of indentation for HTML
        level = dirpath.replace(directory_path, '').count(os.sep)
        indent = '&nbsp;&nbsp;&nbsp;&nbsp;' * level

        # Add current directory as a list item
        html_content += f'<li class="dir">{indent}📁 **{os.path.basename(dirpath) or dirpath}**</li>\n'

        # Add files in the current directory as list items
        for filename in filenames:
            html_content += f'<li class="file">{indent}&nbsp;&nbsp;&nbsp;📄 {filename}</li>\n'

    # End HTML content
    html_content += """
    </ul>
</body>
</html>
"""

    # Write the content to an HTML file
    with open(output_filename, 'w', encoding='utf-8') as f:
        f.write(html_content)

    f.close()
    
    print(f"Successfully generated HTML summary at {os.path.abspath(output_filename)}")

    try:
        linkPath = directory_path + "/" + calendar + ".html"
        relative_src = os.path.relpath( output_filename, os.path.dirname( linkPath))
        print(f"Relative Path: { relative_src }")
        try:
            os.symlink( relative_src, linkPath)
        except:
            pass
    except:
        pass
    # end of eventDaySummary()

def main():
    """
    Main executable for matching transient events
    """
    inDirName = '/home/karl/match/2025/25Jan/25Jan01'
    calendar = '25Jan01'
    outSummaryName = inDirName + "/index.html"
    eventDaySummary('/home/karl/match/2025/25Jan/25Jan01', calendar, outSummaryName)
    
    
if __name__ == "__main__":
    main()

