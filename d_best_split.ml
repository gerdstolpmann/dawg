module Make( L : Loss.LOSS ) = struct

  type best_split_of_features = Feat_map.t -> L.splitter ->
    (Feat.afeature * float * Proto_t.split) option

  let best_split_of_features feature_map splitter =
    Feat_map.fold feature_map (
      fun feature best_opt ->
        let s_opt = splitter#best_split feature in
        match best_opt, s_opt with
          | Some (_, best_loss, best_split), Some (loss, split) ->

            if best_loss < loss then
              (* still superior *)
              best_opt
            else
              (* new champ *)
              Some (feature, loss, split)
          | None, Some (loss, split) ->
            (* first guy's always champ *)
            Some (feature, loss, split)

          | Some _, None -> best_opt
          | None, None -> None

    ) None

end

module LogisticBS = Make(Logistic)
module SquareBS = Make(Square)

type t = [
  | `Logistic of (Logistic.splitter * LogisticBS.best_split_of_features)
  | `Square of (Square.splitter * SquareBS.best_split_of_features)
]