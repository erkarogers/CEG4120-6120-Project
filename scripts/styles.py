# %% [markdown]
# # Shared Styles

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



