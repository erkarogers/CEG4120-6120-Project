# %%
### Import libraries
from ipywidgets import HTML, AppLayout, VBox, Button, Tab, FloatProgress, Output, IntSlider, interactive_output
from ipywidgets import interact, interactive, fixed, interact_manual
import os
import threading
import time
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
    

# %% [markdown]
# ## Create Simulate Tab

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
    value=0.0,
    min=0,
    max=1.0,
    bar_style='info',
    style={'bar_color': '#0000ff'},
    orientation='horizontal'
)

abort_btn = Button(description="Abort", disabled=True)
# abort_btn = ui.RunCommand(label="Run", tooltip="Abort Simulation", start_func=abort)



# https://ipywidgets.readthedocs.io/en/stable/examples/Output%20Widget.html
out = Output(layout={'border': '1px solid black',
                        'width': '100%'})

#@out.capture()
#def function_with_captured_output():
#    print('This goes into the output widget')
    
    
#function_with_captured_output()

#with out:
    #for i in range(5):
     #   print(i, 'Hello world!') 
#        print('what up')
#        print()
    
    

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
header = HTML(
    value="<h2>SEMICONDUCTOR PHYSICS SIMULATION TOOL</h2>",
    placeholder='Title',
    layout=Layout(height='auto')
)
header.add_class('header')

tool_tabs = Tab()
tab_contents = [pw_x_tabs, ph_x_box, epw_x_tabs, simulate_tabs]
tool_tabs.children = tab_contents

tool_tabs.set_title(0, "PW.x Inputs")
tool_tabs.set_title(1, "PH.x Inputs")
tool_tabs.set_title(2, "EPW.x Inputs")
tool_tabs.set_title(3, "Simulate")

app = AppLayout(center=tool_tabs, header=header, pane_heights=[1, 5, '0'])
app.add_class("main")






# %%



