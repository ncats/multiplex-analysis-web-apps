import streamlit as st
import os


def delete_contents_recursively(directory_path):
    for root, dirs, files in os.walk(directory_path, topdown=False):
        for name in files:
            file_path = os.path.join(root, name)
            os.remove(file_path)
            print(f"Deleted file: {file_path}")
        for name in dirs:
            dir_path = os.path.join(root, name)
            os.rmdir(dir_path)
            print(f"Deleted directory: {dir_path}")


def main():

    if st.button("Delete all files in input directory"):
        delete_contents_recursively(os.path.join('.', 'input'))
        st.success("All files deleted.")

    if st.button("Delete all files in output directory"):
        delete_contents_recursively(os.path.join('.', 'output'))
        st.success("All files deleted.")


if __name__ == '__main__':
    main()
