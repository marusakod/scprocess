import re
import os

def define_env(env):
    # --- Existing Software Name Variables ---
    env.variables['sc'] = '<span class="software-name">scprocess</span>'
    env.variables['scsetup'] = '<span class="software-name">scprocess setup</span>'
    env.variables['scnew'] = '<span class="software-name">scprocess newproj</span>'
    env.variables['scrun'] = '<span class="software-name">scprocess run</span>'
    env.variables['scknee'] = '<span class="software-name">scprocess plotknee</span>'

    @env.macro
    def cite_tool(key):
        bib_path = os.path.join(env.project_dir, 'docs/refs.bib')
        
        if not os.path.exists(bib_path):
            return f"Error: refs.bib not found at {bib_path}"

        with open(bib_path, 'r', encoding='utf-8') as f:
            content = f.read()

        escaped_key = re.escape(key)
        pattern = rf"@[a-zA-Z]+\s*\{{{escaped_key}\s*,(.*?)\n\}}"
        match = re.search(pattern, content, re.IGNORECASE | re.DOTALL)

        if not match:
            return f"**Key '{key}' not found**"

        entry_body = match.group(1)

        def get_field(name):
            f_pattern = rf"{name}\s*=\s*[{{|\"|\s]?(.*?)[}}|\"|\s]?,?\n"
            f_match = re.search(f_pattern, entry_body, re.IGNORECASE | re.DOTALL)
            if f_match:
                val = f_match.group(1).strip()
                val = val.replace('{', '').replace('}', '')
                
                if name.lower() == 'author':
                    primary_author = val.split('and')[0].strip()
                    if ',' in primary_author:
                        primary_author = primary_author.split(',')[0]
                    return f"{primary_author} et al."
                return val
            return ""

        author = get_field('author')
        title = get_field('title')
        journal = get_field('journal')
        year = get_field('year')
        note = get_field('note')

        citation = f"{author}. _{title}_. **{journal}**, {year}."
        if note:
            citation += f" {note}"
        
        return citation