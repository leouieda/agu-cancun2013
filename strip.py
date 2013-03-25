"""
Strip the Markdown, newlines and spaces from the abstract in the README
and make it ready for submission (copy pasting).
"""
with open('README.md') as f:
    text = f.readlines()[8:]
text = text[:text.index('## Figures\n')]
print ' '.join(line.strip() for line in text if line != '\n')
