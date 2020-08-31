def isint(var):
    """
    Check if variable is of type `int` or type `float`, but an integer.
    """
    
    if not isinstance(var, int):
        if isinstance(var, float):
            if not var.is_integer():
                return False
            else:
                return True
        else:
            return False
    else:
        return True
            
    