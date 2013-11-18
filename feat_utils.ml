open Dog_t

let array_of_afeature = function
  | `Cat cat -> (
      let categories = Array.of_list cat.c_categories in
      let num_categories = Array.length categories in
      match cat.c_anonymous_category with
        | Some anon_value -> (

            let num_categories = num_categories + 1 in

            let string_opt_of_int value =
              if 0 <= value && value < anon_value then
                Some categories.( value )
              else if value = anon_value then (* anonymous *)
                None
              else if anon_value < value && value < num_categories then
                Some categories.( value - 1 )
              else
                assert false
            in
            match cat.c_vector with
              | `RLE (rle:Rlevec.v) ->
                let result = Array.create rle.Rlevec.length None in
                Rlevec.iter rle (
                  fun ~index ~length ~value ->
                    let string_opt = string_opt_of_int value in
                    for i = index to index + length - 1 do
                      result.(i) <- string_opt
                    done
                );
                `StringAnon result

              | `Dense (vec:Vec.v) ->
                let result = Array.create vec.Vec.length None in
                let width = Utils.num_bytes cat.c_cardinality in
                Vec.iter ~width vec (
                  fun ~index ~value ->
                    let string_opt = string_opt_of_int value in
                    result.(index) <- string_opt
                );
                `StringAnon result
          )
        | None -> (
            let category_0 = List.hd cat.c_categories in
            match cat.c_vector with
              | `RLE rle ->
                let result = Array.create rle.Rlevec.length category_0 in
                Rlevec.iter rle (
                  fun ~index ~length ~value ->
                    let res = categories.( value ) in
                    for i = index to index + length - 1 do
                      result.(i) <- res
                    done

                );
                `String result

              | `Dense vec ->
                let result = Array.create vec.Vec.length category_0 in
                let width = Utils.num_bytes cat.c_cardinality in
                Vec.iter ~width vec (
                  fun ~index ~value ->
                    result.(index) <- categories.( value )
                );
                `String result
          )
    )

  | `Ord { o_vector; o_breakpoints; o_cardinality } -> (
      match o_vector with
        | `RLE rle -> (
            match o_breakpoints with
              | `Float breakpoints ->
                let result = Array.create rle.Rlevec.length 0.0 in
                let breakpoints = Array.of_list breakpoints in

                Rlevec.iter rle (
                  fun ~index ~length ~value ->
                    let res = breakpoints.( value ) in
                    for i = index to index + length - 1 do
                      result.(i) <- res
                    done
                );
                `Float result

              | `Int breakpoints ->
                let result = Array.create rle.Rlevec.length 0 in
                let breakpoints = Array.of_list breakpoints in

                Rlevec.iter rle (
                  fun ~index ~length ~value ->
                    let res = breakpoints.( value ) in
                    for i = index to index + length - 1 do
                      result.(i) <- res
                    done
                );
                `Int result
          )
        | `Dense vec -> (
            let width = Utils.num_bytes o_cardinality in
            match o_breakpoints with
              | `Float breakpoints ->

                let result = Array.create vec.Vec.length 0.0 in
                let breakpoints = Array.of_list breakpoints in
                assert (o_cardinality = Array.length breakpoints);

                Vec.iter ~width vec (
                  fun ~index ~value ->
                    result.( index ) <- breakpoints.( value )
                );
                `Float result

              | `Int breakpoints ->

                let result = Array.create vec.Vec.length 0 in
                let breakpoints = Array.of_list breakpoints in
                assert (o_cardinality = Array.length breakpoints);

                Vec.iter ~width vec (
                  fun ~index ~value ->
                    result.( index ) <- breakpoints.( value )
                );
                `Int result
          )
    )


let id_of_feature = function
  | `Cat { c_feature_id } -> c_feature_id
  | `Ord { o_feature_id } -> o_feature_id

let name_of_feature = function
  | `Cat { c_feature_name_opt } -> c_feature_name_opt
  | `Ord { o_feature_name_opt } -> o_feature_name_opt

let cardinality_of_feature = function
  | `Cat { c_cardinality } -> c_cardinality
  | `Ord { o_cardinality } -> o_cardinality

let vector_of_feature = function
  | `Cat { c_vector } -> c_vector
  | `Ord { o_vector } -> o_vector

let folds_of_feature ~n ~num_folds = function
  | `Ord { o_cardinality; o_vector } ->
    assert ( o_cardinality <= n );
    let cardinality_per_fold = o_cardinality / num_folds in
    if cardinality_per_fold = 0 then
      `TooManyOrdinalFolds o_cardinality
    else
      let folds = Array.create n (-1) in
      (match o_vector with
        | `RLE rle ->
          Rlevec.iter rle (
            fun ~index ~length ~value ->
              let fold = value / cardinality_per_fold in
              for i = index to index + length - 1 do
                folds.(i) <- fold
              done
          );

        | `Dense vec ->
          let width_num_bytes = Utils.num_bytes o_cardinality in
          Vec.iter ~width:width_num_bytes vec (
            fun ~index ~value ->
              let fold = value / cardinality_per_fold in
              folds.(index) <- fold
          );
      );
      `Folds folds

  | `Cat { c_cardinality; c_vector } ->
    if c_cardinality <> num_folds then
      `CategoricalCardinalityMismatch c_cardinality
    else
      let folds = Array.create n (-1) in
      (match c_vector with
        | `RLE rle ->
          Rlevec.iter rle (
            fun ~index ~length ~value ->
              for i = index to index + length - 1 do
                folds.(i) <- value
              done
          );

        | `Dense vec ->
          let width_num_bytes = Utils.num_bytes c_cardinality in
          Vec.iter ~width:width_num_bytes vec (
            fun ~index ~value ->
              folds.(index) <- value
          )
      );
      `Folds folds