import streamlit as st
from streamlit_javascript import st_javascript
from streamlit_extras.add_vertical_space import add_vertical_space 
from st_pages import add_indentation, Page, Section, show_pages

def main():
    st.write('Some Text')
    add_vertical_space(2)

    with st.sidebar:
        url = st_javascript("await fetch('').then(r => window.parent.location.href)")
        st.write(f'''[Open app in new Tab]({url})\n (MS Edge/ Google Chrome)''')

    add_indentation()
    show_pages(
    [
        Page("testingpages.py", "Welcome", "🏠"),
    ])

if __name__ == '__main__':
    main()