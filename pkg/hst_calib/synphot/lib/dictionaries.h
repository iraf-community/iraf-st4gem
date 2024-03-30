# List of compspec commands.  This must be updated if new commands are added
# to the compspec module

string	specdict	"|+|-|x|ebmv|rn|abmag|stmag|vegamag|obmag|mjy\
|jy|fnu|flam|photlam|photnu|counts|bb|pl|hi|grid|"

# List of compband commands.  Must also be updated to include new commands.
string	banddict	"|*|/|=|gauss|lgauss|box|"

# List of forms
string	formdict "|fnu|flam|photnu|photlam|photpix|counts|abmag|stmag|\
vegamag|obmag|jy|mjy|"

# List of other parameters that are considered to be valid "forms" in
# the synphot subroutine
string	lamdict	"|pivlam|rmslam|fwhmlam|barlam|avglam|efflam|"
