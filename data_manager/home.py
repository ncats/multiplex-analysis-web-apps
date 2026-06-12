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
    st.write(f"Welcome {pa.get_current_username()} of {pa.get_user_group(pa.get_current_username())}! 💪")

    st.write("**Note:** If You are encountering issues that seem to be related to large files, memory, or compute power, please try pressing the \"Clear data from memory\" button at left. If you find you constantly need to press this button, please let us know.")

    st.write("We are learning Snowflake's capabilities as you are. If you have trouble using any features, particularly with regard to large files / memory / computer power, please reach out to us and we will solve the problem. Thank you for your patience as we optimize this system!")

    st.space()

    width = 300
    st.page_link(page=st.Page(home.main, url_path='home'), label="Home", icon="🏠", width=width)
    st.page_link(page=st.Page(upload.main, url_path='upload'), label="Upload data", icon="☁️", width=width)
    st.page_link(page=st.Page(download.main, url_path='download'), label="Download data", icon="📥", width=width)
    st.page_link(page=st.Page(delete.main, url_path='delete'), label="Delete data", icon="🔥", width=width)

# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
