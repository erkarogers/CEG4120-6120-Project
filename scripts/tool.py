##
# @mainpage Semiconductor Physics Simulation Tool
# @section description_main Description
# This tool is used to simulate semiconductor physics by using Nanohub and Quantum ESPRESSO.
# @remark The source code documentation was generated using Doxygen.
# @section Table_of_Contents Table of Contents
# - @subpage Tool
# - @subpage EPW_X
# - @subpage PH_X
# - @subpage PW_X
# - @subpage Styles
# 
# @subsection authors Author(s)
# - Source Jupyter Notebook Created by:  Maj. Tim Wolfe
# - GUI Tool Created by:  Ani Bettedpur, Antonio Blackwell, Dylan Clapper, Erin Rogers, Matthew Slusser, Francoise Umubreyi
# - Documentation Last Modified by:  Matthew Slusser on 2023-04-25


# 
# @section notes_main Notes
# 
# @section todo_main TODO
#

##
# @file tool.py
# 
# @brief Main file for running the simulation tool
# 
# @page Tool
# @section description_tool Description
# @ref tool.py is the main driver for the simulation GUI tool
# 
# @section todo_tool TODO
# 
# @section notes_tool Notes
#
# @section libraries_tool Libraries/Modules
# - ipywidgets  : https://ipywidgets.readthedocs.io/en/latest/
#  - HTML       : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#HTML
#  - AppLayout  : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#AppLayout
#  - VBox       : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#VBox
#  - Button     : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Button
#  - Tab        : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Tab
#  - FloatProgress : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#FloatProgress
# - os          : https://docs.python.org/3/library/os.html
# - nanohub     : https://nanohub.org/tools/quantumespresso/
# - epw_x       : epw_x.py
#  - Runs the epw_x.ipynb notebook
# - pw_x        : pw_x.py
#  - Runs the pw_x.ipynb notebook
# - ph_x        : ph_x.py
#  - Runs the ph_x.ipynb notebook
# - styles      : styles.py
#  - Runs the styles.ipynb notebook 
# 


# Import libraries
from ipywidgets import HTML, AppLayout, VBox, Button, Tab, FloatProgress
import os
import threading
import time
%run styles.ipynb
%run pw_x.ipynb 
%run ph_x.ipynb
%run epw_x.ipynb




# %%
# bind Control form to output
## @brief Binds the control form to the output
# @param self The object pointer
# @return None
def on_simulate_btn_click(self):
    bind_PW_X_inputs(self)      # method defined in pw_x.ipynb
    bind_PH_X_inputs(self)      # method defined in ph_x.ipynb
    bind_EPW_X_inputs(self)     # method defined in epw_x.ipynb
    

# %% [markdown]
# ## Create Simulate Tab

# %%
# Simulate
## 
# @brief Simulate button creation
# @param description The description of the button
# @param button_style The style of the button 
# @remark Options: 'success', 'info', 'warning', 'danger' or ''
# @param tooltip The tooltip of the button
# 
simulate_btn = Button(description="Simulate", 
                button_style='', # 'success', 'info', 'warning', 'danger' or '',                    
                tooltip='Run Simulation'
            )
simulate_btn.style.button_color='lightgreen'

simulate_btn.on_click(on_simulate_btn_click)

# simulate_btn = ui.RunCommand(label="Simulate", tooltip="Run Simulation", start_func=start_simulation)

##
# @brief Progress bar creation
# @param value The current value of the progress bar
# @param min The minimum value of the progress bar
# @param max The maximum value of the progress bar
# @param bar_style The style of the progress bar
# @param style Style settings of the progress bar
# @param orientation The orientation of the progress bar
# 
progress_bar = FloatProgress(
    value=0.0,
    min=0,
    max=1.0,
    bar_style='info',
    style={'bar_color': '#0000ff'},
    orientation='horizontal'
)

## 
# @brief Abort button creation
#
abort_btn = Button(description="Abort", disabled=True)
# abort_btn = ui.RunCommand(label="Run", tooltip="Abort Simulation", start_func=abort)

## 
# @brief form_item layout creation
# 
form_items = [
    Box([simulate_btn], layout=form_item_layout()),
    Box([progress_bar], layout=form_item_layout()),
    Box([abort_btn], layout=form_item_layout()),
    Box([out], layout=form_item_layout())
]

simulate_box = Box(form_items, layout=box_layout(40))

#simulate_box


# %% [markdown]
# ## Create Function For Progress Bar Movement

# %%
def bar(progress):
    total = 10
    #progress.value = float(x)/total
    while progressNum<=total:
        progress.value = float(progressNum)/total
        

# %% [markdown]
# ## Create Background Thread for Progress Bar

# %%
#global x
#x=1
#global progressNum
progressNum=0
thread = threading.Thread(target=bar, args=(progress_bar,))

thread.start()

# %% [markdown]
# ## Create Function to Control Status Messages

# %%
statement1='CREATING PW SCF INPUT FILE'
statement2='STARTED PW SCF SIMULATION'
statement3='CREATING PW NSCF INPUT FILE '
statement4='STARTED PW SCF SIMULATION'
statement5='CREATING PH INPUT FILE'
statement6='STARTED PH SIMULATION'
statement7='STARTED POST PROCESSING OF PH SIMULATION RESULTS'
statement8='CREATING EPW INPUT FILE'
statement9='STARTED EPW SIMULATION'
statement10='CREATING EPW OUTPUTS'




global x
x=1
global flag
flag = 0


def printStatements(out):
    total = 10
    #x=0
    while x<=total:
            #print('.')
        if flag == 1:

            if x==1:
                #flag = 0
                t = time.localtime()
                current_time = time.strftime("%H:%M:%S", t)
                with out:
                    print(statement1, ' ',current_time)
                    
    
                reset_flag()  
                
            elif x==2:
                #flag = 0
                t = time.localtime()
                current_time = time.strftime("%H:%M:%S", t)
                with out:
                    print(statement2, ' ',current_time)
                reset_flag()
                
            elif x==3:
                #flag = 0
                t = time.localtime()
                current_time = time.strftime("%H:%M:%S", t)
                with out:
                    print(statement3, ' ',current_time)
                reset_flag()    
                
            elif x==4:
                t = time.localtime()
                current_time = time.strftime("%H:%M:%S", t)
                with out:
                    print(statement4, ' ',current_time)
                reset_flag()
            
            elif x==5:
                t = time.localtime()
                current_time = time.strftime("%H:%M:%S", t)
                with out:
                    print(statement5, ' ',current_time)
                reset_flag()
            
            elif x==6:
                t = time.localtime()
                current_time = time.strftime("%H:%M:%S", t)
                with out:
                    print(statement6, ' ',current_time)
                reset_flag()
            
            elif x==7:
                t = time.localtime()
                current_time = time.strftime("%H:%M:%S", t)
                with out:
                    print(statement7, ' ',current_time)
                reset_flag()
                
            elif x==8:
                t = time.localtime()
                current_time = time.strftime("%H:%M:%S", t)
                with out:
                    print(statement8, ' ',current_time)
                reset_flag()
            
            elif x==9:
                t = time.localtime()
                current_time = time.strftime("%H:%M:%S", t)
                with out:
                    print(statement9, ' ',current_time)
                reset_flag()
            
            elif x==10:
                t = time.localtime()
                current_time = time.strftime("%H:%M:%S", t)
                with out:
                    print(statement10, ' ',current_time)
                reset_flag()
                    
      
def reset_flag():
    global flag
    flag = 0
    #with out:
        #print('changing flag back to 0')

# %% [markdown]
# ## Create Thread for Status Message Updates

# %%
thread2 = threading.Thread(target=printStatements, args=(out,))
thread2.start()

# %% [markdown]
# ## Create Results Tab

# %%
#Results
simulate_btn2 = Button(description="Temporary", 
                button_style='', # 'success', 'info', 'warning', 'danger' or '',                    
                tooltip=''
            )
results_area = Textarea(
    value='Hello World',
    placeholder='Type something',
    description='Output:',
    disabled=True
)
  


form_items = [
    Box([results_area], layout=form_item_layout())
]

results_box = Box(form_items, layout=box_layout(40))





# %% [markdown]
# ## Combine Simulate Tabs

# %%
simulate_tabs = Tab()
tab_contents = [simulate_box,results_box]
simulate_tabs.children = tab_contents

simulate_tabs.set_title(0, "Run")
simulate_tabs.set_title(1, "Results")

# %% [markdown]
# ## Combine All Tabs

# %%
## 
# Combine the all tabs
# @brief This section combines the all tabs from the tool onto one page
# 

## 
# @brief Header of main tool page
# 
header = HTML(
    value="<h2>SEMICONDUCTOR PHYSICS SIMULATION TOOL</h2>",
    placeholder='Title',
    layout=Layout(height='auto')
)
header.add_class('header')

## 
# @brief Tools page tabs creation
tool_tabs = Tab()

## 
# @brief Tools page tabs contents
# @param pw_x_tabs The pw_x page
# @param ph_x_box The ph_x page
# @param epw_x_tabs The epw_x page
# @param simulate_box The simulate page
tab_contents = [pw_x_tabs, ph_x_box, epw_x_tabs, simulate_box]
tool_tabs.children = tab_contents

tool_tabs.set_title(0, "PW.x Inputs")
tool_tabs.set_title(1, "PH.x Inputs")
tool_tabs.set_title(2, "EPW.x Inputs")
tool_tabs.set_title(3, "Simulate")

## 
# @brief Main tool page layout
# 
app = AppLayout(center=tool_tabs, header=header, pane_heights=[1, 5, '0'])
app.add_class("main")






# %%



