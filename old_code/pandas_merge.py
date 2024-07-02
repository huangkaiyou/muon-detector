#%%
import pandas as pd
import numpy as np
array1 =[1,2, 3]
array2 = [5, 6,7]
            
tem = {
        'trash' : [0]
        ,'y' : array1
    }

tem1 = dict(tem)
tem1 = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in tem1.items() ]))
# tem = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in tem.items() ]))

final = {
        
    }

final = dict(final)
final = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in final.items() ]))

# a = pd.concat([final, tem], axis= 1, join='outer')
# %%
