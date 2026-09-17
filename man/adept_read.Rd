\name{adept_read}
\alias{adept_read}
\title{Read an Excel workbook without any external package}
\usage{
adept_read(file_path)
}
\arguments{
\item{file_path}{Path to an \code{.xlsx} (or \code{.csv}) file.}
}
\value{
A named list of data.frames, one per worksheet.
}
\description{
Reads every sheet of an \code{.xlsx} file using only base R. Useful for
inspecting an input file before calling \code{\link{adept}}, and used by the
Shiny front-end for its data preview.
}
\examples{
\dontrun{
sheets <- adept_read("Input.xlsx")
str(sheets[[1]])
}
}
\seealso{
\code{\link{adept}}, \code{\link{adept_gui}}
}
