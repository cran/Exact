updateMat <-
function(Mat, addRow, addedRowVal = 1, nextRowVal = 0) {
  
  for (j in 1:nrow(addRow)) {
    # Add 1 to added row
    Mat[addRow[j,1], addRow[j,2]] <- addedRowVal
    
    # Consider adding 0 to left cell
    if (addRow[j,2] >= 2 && is.na(Mat[addRow[j,1], addRow[j,2]-1])) { Mat[addRow[j,1], addRow[j,2]-1] <- nextRowVal }
    # Consider adding 0 to bottom cell
    if (addRow[j,1] <= (nrow(Mat)-1) && is.na(Mat[addRow[j,1]+1, addRow[j,2]])) { Mat[addRow[j,1]+1, addRow[j,2]] <- nextRowVal }
    # Consider adding 0 to right cell
    if (addRow[j,2] <= (ncol(Mat)-1) && is.na(Mat[addRow[j,1], addRow[j,2]+1])) { Mat[addRow[j,1], addRow[j,2]+1] <- nextRowVal }
    # Consider adding 0 to top cell
    if (addRow[j,1] >= 2 && is.na(Mat[addRow[j,1]-1, addRow[j,2]])) { Mat[addRow[j,1]-1, addRow[j,2]] <- nextRowVal }
  }
  return(Mat)
}
