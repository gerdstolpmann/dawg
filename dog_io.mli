(* write-only *)
module WO : sig
  type t

  val create : string -> int -> t

  val add_feature : t -> Feat.lfeature -> unit

  val close_writer : t -> unit
end

(* read-only *)
module RO : sig
  type t
  val create : string -> t
  val dog : t -> Dog_t.t
  val array : t -> UInt8Array.t

  (* where does the vector sequence end?  This the offset of the last
     byte of the last vector in the sequence, plus one *)
  val end_array_offset : t -> int

end

(* read and append *)
module RW : sig
  type t
  val create : string -> (int * Dog_b.t) option -> t
  (* val create : string -> int -> Dog_t.t -> t *)
  (* [create path size dog_t_] creates a [t], backed by a
     memory-mapped file whose path is [path], and whose size is [size]
     in bytes. *)

  type size_mismatch = {
    expected : int;
    actual : int
  }

  exception SizeMismatch of size_mismatch
  exception FeatureIdNotFound of Dog_t.feature_id

  val write: t -> Dog_t.feature_id -> string -> unit
  (* [write t feature_id vec] writes string [vec] to [t] corresponding
     to [feature_id] . can raise [SizeMismatch] or
     [FeatureIdNotFound]. *)

  val read : t -> Dog_t.feature_id -> string
    (* [read t feature_id] gets the vector corresponding to
       [feature_id], or else raises [FeatureIdNotFound] *)

  val array : t -> UInt8Array.t

  type veq = {
    vector_id : Dog_t.vector_id;
    vector_length : int
  }

  type qfeature = (veq, veq) Dog_t.feature

  val find : t -> Dog_t.feature_id -> qfeature

  val num_observations : t -> int

end
