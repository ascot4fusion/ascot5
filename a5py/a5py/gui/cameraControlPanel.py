'''
Created on Apr 3, 2020

@author: sjjamsa
'''
import tkinter
from .components import NumEntry

def generateCameraControlPanel(tkInterParentElement):
    '''
    Create a tkinter panel for controlling camera initial state.
    @return: dictionary of components in the panel,
             a function that applies the control panel to a a5py.wall.a5vtkwall.camControl -object
    '''
    
    panel  = tkinter.Frame(tkInterParentElement)

    posx   = NumEntry(panel, labeltext="Pos x:", defval=6.0)
    posy   = NumEntry(panel, labeltext="Pos y:", defval=-0.01)
    posz   = NumEntry(panel, labeltext="Pos z:", defval=0.0)

    focx   = NumEntry(panel, labeltext="Foc x:", defval=6.0)
    focy   = NumEntry(panel, labeltext="Foc y:", defval=0.0)
    focz   = NumEntry(panel, labeltext="Foc z:", defval=0.0)

    vux    = NumEntry(panel, labeltext="ViewUp x:", defval=0.0)
    vuy    = NumEntry(panel, labeltext="ViewUp y:", defval=0.0)
    vuz    = NumEntry(panel, labeltext="ViewUp z:", defval=1.0)

    va     = NumEntry(panel, labeltext="VerViewAng:", defval=30.0)

    crMin  = NumEntry(panel, labeltext="ClipRmin:", defval=1.0e-3)
    crMax  = NumEntry(panel, labeltext="ClipRmax:", defval=20.0)


    posx.grid(row=0, column=0, sticky="W")
    posy.grid(row=1, column=0, sticky="W")
    posz.grid(row=2, column=0, sticky="W")

    focx.grid(row=3, column=0, sticky="W")
    focy.grid(row=4, column=0, sticky="W")
    focz.grid(row=5, column=0, sticky="W")

    vux.grid( row=6, column=0, sticky="W")
    vuy.grid( row=7, column=0, sticky="W")
    vuz.grid( row=8, column=0, sticky="W")

    va.grid( row=9, column=0, sticky="W")
    crMin.grid( row=10, column=0, sticky="W")
    crMax.grid( row=11, column=0, sticky="W")

    components =  {'panel':panel,
            'posx':posx, 'posy':posy, 'posz':posz,
            'focx':focx, 'focy':focy, 'focz':focz,
             'vux': vux,  'vuy': vuy,  'vuz': vuz,
             'va' : va , 'crMin':crMin, 'crMax':crMax             }
    
    def applyCameraControlPanel(components,camControl):

        #import a5py.wall
        #camControl=a5py.wall.a5vtkwall()
        
        camControl.position   = ( components['posx'].getval(), 
                                  components['posy'].getval(), 
                                  components['posz'].getval()  )
        
        camControl.focalPoint = ( components['focx'].getval(), 
                                  components['focy'].getval(), 
                                  components['focz'].getval()  )
        
        camControl.viewUp     = ( components[ 'vux'].getval(), 
                                  components[ 'vuy'].getval(), 
                                  components[ 'vuz'].getval()  )

        camControl.viewAngle  =   components[ 'va' ].getval() 

        camControl.clipRange = ( components['crMin'].getval(),
                                 components['crMax'].getval()   )

    return components, applyCameraControlPanel
    
