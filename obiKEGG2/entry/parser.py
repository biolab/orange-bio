"""
A parser for DBGET entries
 
"""

class DBGETEntryParser(object):
    """ A DBGET entry parser (inspired by ``xml.dom.pulldom``)
    """
    ENTRY_START = 0
    ENTRY_END = 1
    SECTION_START = 2
    SECTION_END = 3
    SUBSECTION_START = 4
    SUBSECTION_END = 5
    TEXT = 6
    
    def __init__(self):
        pass
    
    def parse(self, stream):
        entry_offset = None
        section_offset = None
        section_title = None
        subsection_offset = None
        subsection_title = None
        textline_start = None
        
        for line in stream:
            startswith = line.startswith
            # TODO: Reorder by frequency 
            if startswith("ENTRY"):
                # Start parsing new entry
                yield (self.ENTRY_START, None, None)
                title, rest = self._partition_section_title(line)
                entry_offset = len(line) - len(rest)
                textline_start = " " * entry_offset
                yield (self.SECTION_START, title, rest)
                yield (self.SECTION_END, title, None)
                
            elif startswith("///"):
                # End entry
                if subsection_title is not None:
                    # End current subsection if any
                    yield (self.SUBSECTION_END, subsection_title, None)
                    subsection_title = None
                    subsection_offset = None
                    
                if section_title is not None:
                    # End current section if any
                    yield (self.SECTION_END, section_title, None)
                    section_title = None
                    subsection_offset = None
                    
                yield (self.ENTRY_END, None, None)
                entry_offset = None
                textline_start = None
                
            elif not startswith(" "):
                # Start new section
                if subsection_title is not None:
                    # End current subsection if any
                    yield (self.SUBSECTION_END, subsection_title, None)
                    subsection_title = None
                    subsection_offset = None
                    
                if section_title is not None:
                    # End current section if any
                    yield (self.SECTION_END, section_title, None)
                    section_title = None
                    subsection_offset = None
                title, rest = self._partition_section_title(line)
                section_offset = len(line) - len(rest)
                section_title = title
                yield (self.SECTION_START, section_title, rest)
                
            elif startswith(textline_start):
                # A line of text
                # TODO: pass the current subsection/section title
                yield (self.TEXT, None, line[entry_offset:])
                
            elif startswith(" "):
                # Start a new subsection
                if subsection_title is not None:
                    # End current subsection
                    yield (self.SUBSECTION_END, subsection_title, None)
                title, rest = self._partition_subsection_title(line)
                subsection_offset = len(line) - len(rest)
                subsection_title = title
                yield (self.SUBSECTION_START, subsection_title, rest)
                
        # Close any remaining sections/entries
        if subsection_title is not None:
            yield (self.SUBSECTION_END, subsection_title, None)
        if section_title is not None:
            yield (self.SECTION_END, section_title, None)
        if entry_offset is not None:
            yield (self.ENTRY_END, None, None)
    
    def parse_string(self, string):
        from StringIO import StringIO
        return self.parse(StringIO(string))
    
    def _partition_section_title(self, line):
        """ Split the section title from the rest of the line
        """
        title, rest = line.split(" ", 1)
        rest = rest.lstrip(" ")
        return title, rest
    
    def _partition_subsection_title(self, line):
        """ Split the subsection title from the rest of the line
        """
        line = line.lstrip(" ")
        return self._partition_section_title(line)
    