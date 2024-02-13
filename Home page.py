import streamlit as st

from PIL import Image


from plotly.subplots import make_subplots
import plotly.graph_objects as go


st.title('Welcome to Virialis!')


st.write("Below, you can see a diagram on how to use Virialis")
image = Image.open('pics/1.png')

st.image(image)