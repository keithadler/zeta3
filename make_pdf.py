"""Generate paper.pdf from paper.md with embedded figures and math rendering."""

from weasyprint import HTML
import markdown
import re

with open('paper.md', 'r') as f:
    md_content = f.read()

# Convert display math $$ ... $$ to HTML before markdown processing
# Replace $$ ... $$ with a styled div
def replace_display_math(match):
    expr = match.group(1).strip()
    return f'<div class="math-display">{expr}</div>'

def replace_inline_math(match):
    expr = match.group(1).strip()
    return f'<span class="math-inline">{expr}</span>'

# Handle display math first ($$...$$)
md_content = re.sub(r'\$\$(.*?)\$\$', replace_display_math, md_content, flags=re.DOTALL)

# Handle inline math ($...$) - be careful not to match currency
md_content = re.sub(r'(?<!\$)\$([^\$\n]+?)\$(?!\$)', replace_inline_math, md_content)

html_content = markdown.markdown(md_content, extensions=['tables', 'fenced_code'])

# Add MathJax-like rendering via CSS (since we can't use JS in PDF)
# We'll use a monospace italic font for math expressions
full_html = '''<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
body { font-family: 'Times New Roman', serif; max-width: 800px; margin: 40px auto; padding: 0 20px; line-height: 1.6; font-size: 11pt; }
h1 { font-size: 18pt; text-align: center; margin-bottom: 5px; }
h2 { font-size: 14pt; border-bottom: 1px solid #ccc; padding-bottom: 5px; margin-top: 30px; page-break-before: auto; }
h3 { font-size: 12pt; margin-top: 20px; }
h4 { font-size: 11pt; }
table { border-collapse: collapse; width: 100%; margin: 15px 0; font-size: 9pt; }
th, td { border: 1px solid #ddd; padding: 6px 8px; text-align: left; }
th { background-color: #f5f5f5; font-weight: bold; }
code { background-color: #f4f4f4; padding: 2px 4px; font-size: 9pt; }
pre { background-color: #f4f4f4; padding: 12px; overflow-x: auto; font-size: 9pt; }
blockquote { border-left: 3px solid #333; padding-left: 15px; margin-left: 0; font-weight: bold; }
img { max-width: 85%; height: auto; margin: 15px auto; display: block; page-break-inside: avoid; }
em { font-style: italic; }
strong { font-weight: bold; }
hr { border: none; border-top: 1px solid #ccc; margin: 30px 0; }
.math-display { font-family: 'Times New Roman', serif; font-style: italic; text-align: center; margin: 20px 0; padding: 10px; font-size: 11pt; }
.math-inline { font-family: 'Times New Roman', serif; font-style: italic; }
</style>
</head>
<body>
''' + html_content + '''
</body>
</html>'''

HTML(string=full_html, base_url='.').write_pdf('paper.pdf')
print('Done: paper.pdf')
