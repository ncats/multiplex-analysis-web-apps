# Import relevant libraries.
import streamlit as st
import upload
import download
import delete
import home
import framework.platform_abstraction as pa


# Define the main function.
def main():

    # Welcome message.
    st.write(f"Welcome {pa.get_current_username()} of {pa.get_user_group(pa.get_current_username())}! 💪 Please navigate below:")

    st.space()

    width = 300
    st.page_link(page=st.Page(home.main, url_path='home'), label="Home", icon="🏠", width=width)
    st.page_link(page=st.Page(upload.main, url_path='upload'), label="Upload data", icon="☁️", width=width)
    st.page_link(page=st.Page(download.main, url_path='download'), label="Download data", icon="📥", width=width)
    st.page_link(page=st.Page(delete.main, url_path='delete'), label="Delete data", icon="🔥", width=width)

# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
