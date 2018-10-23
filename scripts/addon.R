exclusive_or <- function(a, b){
  return(setdiff(union(a,b), intersect(a,b)))
}

exclusive_or_3 <- function(a,b,c){
  abc = union(a,b)
  abc = union(abc,c)
  
  ab = intersect(a,b)
  bc = intersect(b,c)
  ac = intersect(a,c)
  re = setdiff(abc, ab)
  re = setdiff(re, bc)
  re = setdiff(re, ac)
  return(re)
}