(add-hook 'c++-mode-hook
  (function
    (lambda ()
  ;; SCIP customizations for c-mode and c++-mode
  (setq-default c-basic-offset 3)
  (c-set-offset 'substatement-open 0)
  (c-set-offset 'statement-case-open 0)
  (c-set-offset 'brace-list-open '-)
  (c-set-offset 'inextern-lang '0)
  (c-set-offset 'arglist-intro '+)
  (c-set-offset 'arglist-cont 0)
  (c-set-offset 'arglist-cont-nonempty '+)
  (c-set-offset 'arglist-close '+)
  (set-variable 'fill-column 120)
 ;; this will make sure spaces are used instead of tabs
  (setq tab-width 8 indent-tabs-mode nil)
  )))
