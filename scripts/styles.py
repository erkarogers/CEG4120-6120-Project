# %% [markdown]
# # Shared Styles

##
# @file styles.py
# 
# @brief Shared styles for the GUI tool
# 
# @page Styles
# @section description_styles Description
# @ref styles.py contains the shared styles for the GUI tool
# 
# @section todo_styles TODO
# 
# @section notes_styles Notes
#
# @section libraries_styles Libraries/Modules
# - ipywidgets  : https://ipywidgets.readthedocs.io/en/latest/
#  - Layout     : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Layout
# 


# %%
from ipywidgets import Layout

# %%
def input_layout(width):
    return Layout(width='{}%'.format(width))

# %%
def form_item_layout():
    return Layout(
    display='flex',
    flex_flow='row',
    justify_content='space-between')
    #justify_content='center')
    #justify_content='space-around')

# %%
def box_layout(width):
    return Layout(
    display='flex',
    flex_flow='column',
    #border='solid 2px',
    align_items='stretch',
    width='{}%'.format(width)
)

# %%
%%html
<style>
.header{
    background-color:#dedede;
    padding: 10px;
}
.main {
    border: solid 1px;
}

.tab {
    border: 0px;
}
</style>

# %%



