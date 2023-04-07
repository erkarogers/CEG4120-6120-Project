# %%
### Import libraries
from ipywidgets import HTML, AppLayout, VBox, Button, Tab, FloatProgress
import os
%run styles.ipynb
%run pw_x.ipynb 
%run ph_x.ipynb
%run epw_x.ipynb

# %%
# bind Control form to output

def on_simulate_btn_click(self):
    bind_PW_X_inputs(self)      # method defined in pw_x.ipynb
    bind_PH_X_inputs(self)      # method defined in ph_x.ipynb
    bind_EPW_X_inputs(self)     # method defined in epw_x.ipynb
    

# %%
# Simulate
simulate_btn = Button(description="Simulate", 
                button_style='', # 'success', 'info', 'warning', 'danger' or '',                    
                tooltip='Run Simulation'
            )
simulate_btn.style.button_color='lightgreen'
simulate_btn.on_click(on_simulate_btn_click)

# simulate_btn = ui.RunCommand(label="Simulate", tooltip="Run Simulation", start_func=start_simulation)
progress_bar = FloatProgress(
    value=4.5,
    min=0,
    max=10.0,
    bar_style='info',
    style={'bar_color': '#0000ff'},
    orientation='horizontal'
)

abort_btn = Button(description="Abort", disabled=True)
# abort_btn = ui.RunCommand(label="Run", tooltip="Abort Simulation", start_func=abort)

form_items = [
    Box([simulate_btn], layout=form_item_layout()),
    Box([progress_bar], layout=form_item_layout()),
    Box([abort_btn], layout=form_item_layout())
]

simulate_box = Box(form_items, layout=box_layout(40))

#simulate_box


# %% [markdown]
# ## Combine All Tabs

# %%
header = HTML(
    value="<h2>SEMICONDUCTOR PHYSICS SIMULATION TOOL</h2>",
    placeholder='Title',
    layout=Layout(height='auto')
)
header.add_class('header')

tool_tabs = Tab()
tab_contents = [pw_x_tabs, ph_x_box, epw_x_tabs, simulate_box]
tool_tabs.children = tab_contents

tool_tabs.set_title(0, "PW.x Inputs")
tool_tabs.set_title(1, "PH.x Inputs")
tool_tabs.set_title(2, "EPW.x Inputs")
tool_tabs.set_title(3, "Simulate")

app = AppLayout(center=tool_tabs, header=header, pane_heights=[1, 5, '0'])
app.add_class("main")


# %%



