(* GS. Finding bins by analyzing frequencies. The algorithm first creates
   b bins in a trivial way, then adds more elements one by one, and in every
   step two adjacent bins are merged, so that we have before the step and 
   after the step exactly b bins.

   Which bins are merged? The two bins with the lowest combined counts.
   This is analogous to Huffman coding, where in one step the nodes with the
   lowest weights are combined into a bit.

   By design, this works only well if there are initial counts. If all counts are
   just 1 (which can easily happen for inputs coming from real measurements)
   the algorithm just pairs up k adjacent elements into a single one. (I guess
   in that case we'd better off with a compressor that also takes the value
   differences into account - bins closely together are more likely to paired
   up than those lying far away from each other.)
 *)


type bin = {
  left: float;
  mean: float;
  right : float
}

let rec huf_loop best best_f prefix = function
  | [] | [_] ->
    let (best_prefix, best_node, best_tail) = best in
    List.rev_append best_prefix (best_node :: best_tail)
  | hd :: ((nxt :: nxt_tl) as tl) ->
    let (hd_bin, hd_f) = hd in
    let (nxt_bin, nxt_f) = nxt in
    let candidate_f = hd_f + nxt_f in
    if candidate_f >= best_f then
      huf_loop best best_f (hd :: prefix) tl
    else
      let total_f = float (hd_f + nxt_f) in
      let mean =
        (hd_bin.mean *. float hd_f +. nxt_bin.mean *. float nxt_f) /. total_f
      in
      let candidate_bin = {
        left = hd_bin.left;
        mean;
        right = nxt_bin.right;
      } in
      let candidate_node = candidate_bin, candidate_f in
      let candidate = (prefix, candidate_node, nxt_tl) in
      huf_loop candidate candidate_f (hd :: prefix) tl

let huf_one b_plus_one =
  match b_plus_one with
    | hd :: ((nxt :: nxt_tl) as tl) ->
      let (hd_bin, hd_f) = hd in
      let (nxt_bin, nxt_f) = nxt in
      let candidate_f = hd_f + nxt_f in
      let total_f = float (hd_f + nxt_f) in
      let mean =
        (hd_bin.mean *. float hd_f +. nxt_bin.mean *. float nxt_f) /. total_f
      in
      let candidate_bin = {
        left = hd_bin.left;
        mean;
        right = nxt_bin.right;
      } in
      let candidate_node = candidate_bin, candidate_f in
      let candidate = ([], candidate_node, nxt_tl) in
      huf_loop candidate candidate_f [hd] tl
    | _ -> b_plus_one

let singleton x =
  { left = x; mean = x; right = x }

(* GS. Better done with two passes, one extracting the n first elements, and
   one skipping the n first elements. Doesn't stress the GC as much.
 *)
let list_split_n n l =
  let rec iter n prefix tail = match tail with
    | _ when n <= 0 -> prefix, tail
    | [] -> prefix, tail
    | hd :: tl -> iter (pred n) (hd :: prefix) tl
  in
  iter n [] l

(* GS:
    b = number of bins
    elem_counts = list of (value,count), sorted by descending values
 *)

let create elem_counts b =
  let b_elems, more_elems = list_split_n b elem_counts in
  let init_hist = List.map (fun (x, f) -> singleton x, f) b_elems in
  List.fold_left (fun hist (x, f) ->
    let b_plus_one = (singleton x, f) :: hist in
    huf_one b_plus_one
  ) init_hist more_elems
