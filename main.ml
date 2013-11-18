open Cmdliner

let _ =
  let main_cmd =
    let doc = "the EigenDog stochastic gradient boosted tree learner" in
    Term.(pure ()), Term.info "dawg" ~doc
  in

  let commands = List.flatten [
      Csv2dog.commands;
      Dog_inspect.commands;
      Learn.commands;
      Eval.commands;
      Bi_metrics.commands;
      Gen_code.commands;
    ] in

  match Term.eval_choice ~catch:false main_cmd commands with
    | `Error _ -> exit 1
    | _ -> exit 0