# .lintr.R
linters <- linters_with_defaults(
  line_length_linter        = NULL,
  object_name_linter        = NULL,
  commented_code_linter     = NULL,
  object_length_linter      = NULL,
  return_linter             = NULL,
  assignment_linter         = NULL,
  comma_linter              = NULL,
  cyclocomp_linter          = NULL,
  duplicate_linter          = NULL,
  function_name_linter      = NULL,
  infix_spaces_linter       = NULL,
  no_tab_linter             = NULL,
  object_usage_linter       = NULL,
  semicolon_terminator_linter = NULL,
  single_quotes_linter      = NULL,
  spaces_inside_linter      = NULL,
  T_and_F_symbol_linter     = NULL,
  trailing_blank_lines_linter = NULL,
  trailing_whitespace_linter  = NULL
)