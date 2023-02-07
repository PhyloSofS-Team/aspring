import os

def get_immediate_subdirectories(a_dir): #while inside directory get subdirectories names into list
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]
