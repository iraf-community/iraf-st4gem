# Define the Font linked list memory management.
define  FONT_ENTRY_SIZE 2
define  FONT_PREVIOUS   Memi[$1]
define  FONT_ENTRY_PTR  Memi[$1+1]
define  FONT_ENTRY      Memc[FONT_ENTRY_PTR($1)]

#---------------------------------------------------------------------------
.help psk_addfont 1May92 plot
.ih
NAME
psk_addfont -- Add a font name to a font list.
psk_popfont -- Retrieve/remove a font from the list.
.ih
USAGE
call psk_addfont (font, font_list)
call psk_popfont (font_list, font, max_size)
.ih
ARGUMENTS
.ls font (char[ARB])
The font name to add to the list.
.le
.ls font_list (pointer)
Pointer to first entry in the font list.
.le
.ls max_size (int)
The maximum length of a font name to return.
.le
.ih
DESCRIPTION
Since fonts can be added arbitrarily during drawing, we need to keep a list
of changes to add to the PostScript at the end of the file.  All duplications
are ignored.
.endhelp
#---------------------------------------------------------------------------

procedure psk_addfont (font, font_list)

char    font[ARB]                  # I:  Name to add to list.
pointer font_list                  # IO: The font list.

# Declarations.
bool    found                      # True if font is already in list.

pointer entry                      # Pointer to font entries.

# Function prototypes.
bool    streq()

begin
        
        # Search the list to make sure that the font hasn't already
        # been specified.
        entry = font_list
        found = false
        while (entry != NULL && !found) {
            found = streq (FONT_ENTRY(entry), font)
            entry = FONT_PREVIOUS(entry)
        }
        
        # If it isn't found, add the font to the list.
        if (!found) {
            
            entry = font_list
            call calloc (font_list, FONT_ENTRY_SIZE, TY_STRUCT)
            
            FONT_PREVIOUS(font_list) = entry
            call calloc(FONT_ENTRY_PTR(font_list), SZ_LINE, TY_CHAR)
            call strcpy (font, FONT_ENTRY(font_list), SZ_LINE)
            
        }
        
end
#---------------------------------------------------------------------------
# End of psk_addfont
#---------------------------------------------------------------------------
# psk_popfont -- Retrieve and remove a font from the font list.
#---------------------------------------------------------------------------

procedure psk_popfont (font_list, font, max_size)

pointer font_list                  # IO: The font list.
char    font[max_size]             # O:  The next font off the list.
int     max_size                   # I:  The maximum font name size.

# Declarations
pointer current_entry              # The current font entry.

begin
        
        if (font_list == NULL)
            font[1] = EOS
        else {
            current_entry = font_list
            call strcpy (FONT_ENTRY(current_entry), font, max_size)
            call mfree (FONT_ENTRY_PTR(current_entry), TY_CHAR)
            font_list = FONT_PREVIOUS(current_entry)
            call mfree (current_entry, TY_STRUCT)
        }
        
end
#---------------------------------------------------------------------------
# End of psk_popfont
#---------------------------------------------------------------------------
