"""
https://discourse.julialang.org/t/enums-to-string-string-to-enum-a-simple-guide-to-use-enums-strings/126044
"""
string_to_enum(str::String) = eval(Symbol(str))

"""
https://discourse.julialang.org/t/enums-to-string-string-to-enum-a-simple-guide-to-use-enums-strings/126044
"""
enum_to_string(enum_value::Enum) = string(Symbol(enum_value))
