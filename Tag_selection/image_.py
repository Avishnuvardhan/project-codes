from numpy.core.fromnumeric import size
import streamlit as st
import pandas as pd
import numpy as np
from PIL import Image
import os 

def get_data(csv_path):
    df = pd.read_csv(csv_path)
    df = df.replace(np.nan, '')

    # Select Color of the Tag
    color = df['Color'].drop_duplicates()
    Select_Color = st.sidebar.selectbox('Select the Color of Tag',color)

    # Select Shape of the Tag
    shape = df['Shape'].loc[df['Color']==Select_Color]
    Select_Shape = st.sidebar.selectbox("Select Shape of the tag",shape)

    # Select Tag Size
    size = df['Tag_Size'].loc[(df['Color']==Select_Color)& (df['Shape']==Select_Shape)]
    Select_Tag = st.sidebar.selectbox("Select Tag Size", size)

    # df.loc[(df['Color']==color) & (df['Shape']==shape) & (df['Tag_Size']==Tag)]
    Image_url = df.loc[(df['Color']==Select_Color)& (df['Shape']==Select_Shape) ]
    # st.write("Results",Image_url )

    Images_names = []
    # Display images
    for image_name in list(Image_url['Tag_Name']):
        storage_path = r'C:\Users\user\Desktop\me\Data'
        storage_path = storage_path + r"\{0}".format(image_name)
        Images_names.append(storage_path)

    return Images_names

if __name__ == '__main__':
     Images_names = get_data(r'C:\Users\user\Desktop\me\new_excel.csv')
     
     for k in Images_names:
        img = Image.open(k)
        if img.mode != 'RGB':
            img = img.convert('RGB')
            st.image(img,width=100)
        else:
            st.image(img,width=100)