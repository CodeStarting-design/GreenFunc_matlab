import PyPDF2
reader = PyPDF2.PdfReader('s11432-009-0021-0.pdf')
text = ""
for page in reader.pages:
    text += page.extract_text()
with open('pdf_text.txt', 'w') as f:
    f.write(text)
