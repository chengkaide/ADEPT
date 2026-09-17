# ADEPT: Automated Depth Profiling Technique
# xlsx.R — minimal .xlsx reader and writer in pure base R
#
# Removes the readxl / openxlsx / writexl dependencies. An .xlsx file is a ZIP
# archive holding OOXML parts; this module reads those parts with `unzip()` and
# regex parsing, and writes them back as a ZIP with *stored* (uncompressed)
# entries, which needs no deflate library.
#
# Reading is validated against readxl on the bundled example workbook; writing
# is validated by reading the output back with readxl.

# ---------------------------------------------------------------------------
#  Small helpers
# ---------------------------------------------------------------------------

#' Convert Excel column letters (A, B, ..., AA) to a 1-based index
#'
#' Lookup-table based: building the reference table once turns the previous
#' per-cell `vapply` + `utf8ToInt` into a single `match()` call, which matters
#' on workbooks with >100k cells.
#' @keywords internal
ADEPT_COL_REFS <- local({
  n <- 18278L                      # 1 .. ZZZ
  out <- character(n)
  for (i in seq_len(n)) {
    v <- i; s <- ""
    while (v > 0L) {
      r <- (v - 1L) %% 26L
      s <- paste0(intToUtf8(65L + r), s)
      v <- (v - 1L) %/% 26L
    }
    out[i] <- s
  }
  out
})

#' @keywords internal
col_letters_to_int <- function(x) {
  match(x, ADEPT_COL_REFS)
}

#' Escape text for XML
#' @keywords internal
xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&apos;", x, fixed = TRUE)
  x
}

#' Unescape XML text, including numeric character references
#' @keywords internal
xml_unescape <- function(x) {
  x <- gsub("&lt;", "<", x, fixed = TRUE)
  x <- gsub("&gt;", ">", x, fixed = TRUE)
  x <- gsub("&quot;", "\"", x, fixed = TRUE)
  x <- gsub("&apos;", "'", x, fixed = TRUE)
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  x <- gsub("_x000D_", "\r", x, fixed = TRUE)
  x
}

#' Read a whole file as one UTF-8 string
#' @keywords internal
read_utf8 <- function(path) {
  n <- file.info(path)$size
  if (is.na(n) || n == 0) return("")
  con <- file(path, open = "rb", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  raw <- readBin(con, "raw", n)
  if (length(raw) >= 3 && raw[1] == as.raw(0xEF) && raw[2] == as.raw(0xBB) &&
      raw[3] == as.raw(0xBF)) raw <- raw[-(1:3)]
  enc2utf8(rawToChar(raw))
}

#' Read one member of a ZIP archive directly, without extracting to disk
#'
#' Uses `unz()` so no temporary files are created (and none have to be
#' deleted, which is unreliable on Windows).
#'
#' @param zipfile Path to the archive.
#' @param entry Member name, e.g. `"xl/worksheets/sheet1.xml"`.
#' @return The member content as one UTF-8 string, or `""` if absent.
#' @keywords internal
read_zip_entry <- function(zipfile, entry) {
  con <- try(unz(zipfile, entry, open = "rb", encoding = "UTF-8"),
             silent = TRUE)
  if (inherits(con, "try-error")) return("")
  on.exit(close(con), add = TRUE)
  raw <- try(readBin(con, "raw", 1e9), silent = TRUE)
  if (inherits(raw, "try-error") || length(raw) == 0) return("")
  if (length(raw) >= 3 && raw[1] == as.raw(0xEF) && raw[2] == as.raw(0xBB) &&
      raw[3] == as.raw(0xBF)) raw <- raw[-(1:3)]
  enc2utf8(rawToChar(raw))
}

#' Names of all members of a ZIP archive
#' @keywords internal
zip_entries <- function(zipfile) {
  info <- try(utils::unzip(zipfile, list = TRUE), silent = TRUE)
  if (inherits(info, "try-error")) return(character(0))
  as.character(info$Name)
}

# ---------------------------------------------------------------------------
#  Reading
# ---------------------------------------------------------------------------

#' Extract the text of every <si> entry in sharedStrings.xml
#' @keywords internal
parse_shared_strings <- function(xml) {
  if (!nzchar(xml)) return(character(0))
  items <- strsplit(xml, "<si>", fixed = TRUE)[[1]][-1]
  if (length(items) == 0) return(character(0))
  items <- sub("</si>.*$", "", items)
  m <- gregexpr("<t[^>]*>(.*?)</t>", items, perl = TRUE)
  vapply(regmatches(items, m), function(parts) {
    if (length(parts) == 0) return("")
    parts <- sub("^<t[^>]*>", "", parts)
    parts <- sub("</t>$", "", parts)
    xml_unescape(paste(parts, collapse = ""))
  }, character(1))
}

#' Parse one worksheet part into a character matrix
#'
#' @param xml Worksheet XML text.
#' @param shared Character vector of shared strings.
#' @return List with `nrow`, `ncol` and `cells` (column-major character vector).
#' @keywords internal
parse_worksheet <- function(xml, shared) {
  body <- sub("^.*?<sheetData[^>]*>", "", xml)
  body <- sub("</sheetData>.*$", "", body)

  chunks <- strsplit(body, "<c ", fixed = TRUE)[[1]][-1]
  if (length(chunks) == 0) {
    return(list(nrow = 0L, ncol = 0L, cells = character(0)))
  }

  ref  <- sub('^r="([A-Z]+)([0-9]+)".*$', "\\1", chunks)
  rown <- suppressWarnings(as.integer(sub('^r="[A-Z]+([0-9]+)".*$', "\\1", chunks)))
  # Cells without an explicit reference are rare; fall back to sequential fill.
  if (all(is.na(rown))) {
    rown <- seq_along(chunks)
    ref  <- rep(NA_character_, length(chunks))
  }
  coln <- col_letters_to_int(ref)
  if (anyNA(coln)) {
    miss <- which(is.na(coln))
    for (i in miss) coln[i] <- if (i > 1 && !is.na(coln[i - 1])) coln[i - 1] + 1L else 1L
  }

  is_shared   <- grepl('t="s"', chunks, fixed = TRUE)
  is_inline   <- grepl('t="inlineStr"', chunks, fixed = TRUE)
  is_str      <- grepl('t="str"', chunks, fixed = TRUE)
  has_v       <- grepl("<v>", chunks, fixed = TRUE)

  val <- rep(NA_character_, length(chunks))

  if (any(has_v)) {
    vv <- sub("^.*?<v>(.*?)</v>.*$", "\\1", chunks[has_v], perl = TRUE)
    idx <- which(has_v)
    if (any(is_shared[idx])) {
      # sharedStrings indices are 0-based; R vectors are 1-based.
      si <- suppressWarnings(as.integer(vv)) + 1L
      sel <- is_shared[idx]
      vals <- shared[si[sel]]
      vals[is.na(si[sel]) | si[sel] < 1L | si[sel] > length(shared)] <- NA
      vv[sel] <- vals
    }
    if (any(is_str[idx])) vv[is_str[idx]] <- xml_unescape(vv[is_str[idx]])
    if (any(is_inline[idx])) vv[is_inline[idx]] <- xml_unescape(vv[is_inline[idx]])
    val[idx] <- vv
  }

  if (any(is_inline & !has_v)) {
    ii <- which(is_inline & !has_v)
    val[ii] <- xml_unescape(sub("^.*?<is>(.*?)</is>.*$", "\\1",
                                chunks[ii], perl = TRUE))
    val[ii] <- gsub("</?t[^>]*>", "", val[ii])
  }

  nr <- if (length(rown)) max(rown, na.rm = TRUE) else 0L
  nc <- if (length(coln)) max(coln, na.rm = TRUE) else 0L
  cells <- rep(NA_character_, nr * nc)
  cells[(coln - 1L) * nr + rown] <- val

  list(nrow = nr, ncol = nc, cells = cells)
}

#' Turn a parsed worksheet into a data.frame with header + type guessing
#' @keywords internal
worksheet_to_df <- function(ws) {
  if (ws$nrow < 1 || ws$ncol < 1) {
    return(data.frame())
  }
  m <- matrix(ws$cells, nrow = ws$nrow, ncol = ws$ncol)
  header <- as.character(m[1, ])
  header[is.na(header) | !nzchar(header)] <-
    paste0("V", which(is.na(header) | !nzchar(header)))
  m <- m[-1, , drop = FALSE]

  out <- vector("list", ncol(m))
  for (j in seq_len(ncol(m))) {
    col <- m[, j]
    num <- suppressWarnings(as.numeric(col))
    if (sum(!is.na(num)) >= sum(!is.na(col)) && any(!is.na(num))) {
      out[[j]] <- num
    } else {
      out[[j]] <- col
    }
  }
  names(out) <- make.unique(header)
  n <- if (nrow(m)) nrow(m) else 0L
  as.data.frame(lapply(out, function(v) v[seq_len(n)]),
                stringsAsFactors = FALSE, optional = TRUE)
}

#' Read all sheets of an .xlsx file without any external package
#'
#' @param path Path to the workbook.
#' @return Named list of data.frames, one per sheet.
#' @keywords internal
#' Read an Excel workbook without any external package
#'
#' Reads every sheet of an \code{.xlsx} file using only base R. Useful for
#' inspecting an input file before calling \code{\link{adept}}, and used by the
#' Shiny front-end for its data preview.
#'
#' @param file_path Path to an \code{.xlsx} (or \code{.csv}) file.
#' @return A named list of data.frames, one per worksheet.
#'
#' @seealso \code{\link{adept}}, \code{\link{adept_gui}}
#' @export
#' @examples
#' \dontrun{
#' sheets <- adept_read("Input.xlsx")
#' str(sheets[[1]])
#' }
adept_read <- function(file_path) {
  read_input(file_path)
}

#' Read all sheets of an .xlsx file without any external package
#'
#' @param path Path to the workbook.
#' @return Named list of data.frames, one per sheet.
#' @keywords internal
read_xlsx_base <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)

  entries <- zip_entries(path)
  if (length(entries) == 0) {
    stop("Could not list the contents of '", basename(path),
         "'. Is it a valid .xlsx (ZIP) file?", call. = FALSE)
  }

  pick <- function(want) {
    hit <- entries[tolower(entries) == tolower(want)]
    if (length(hit) == 0) "" else hit[1]
  }

  read_part <- function(want) {
    nm <- pick(want)
    if (!nzchar(nm)) return("")
    read_zip_entry(path, nm)
  }

  wb_xml   <- read_part("xl/workbook.xml")
  rels_xml <- read_part("xl/_rels/workbook.xml.rels")
  ss_xml   <- read_part("xl/sharedStrings.xml")

  shared <- parse_shared_strings(ss_xml)

  # sheet name -> target part
  sheets <- regmatches(
    wb_xml,
    gregexpr('<sheet[[:space:]][^>]*>', wb_xml, perl = TRUE)
  )[[1]]
  sheet_names <- sub('^.*?name="([^"]*)".*$', "\\1", sheets)
  sheet_rids  <- sub('^.*?r:id="([^"]*)".*$', "\\1", sheets)

  rels <- regmatches(
    rels_xml,
    gregexpr('<Relationship[[:space:]][^>]*>', rels_xml, perl = TRUE)
  )[[1]]
  rel_ids  <- sub('^.*?Id="([^"]*)".*$', "\\1", rels)
  rel_targ <- sub('^.*?Target="([^"]*)".*$', "\\1", rels)

  if (length(sheet_names) == 0) {
    # No workbook.xml (unusual) — fall back to the first worksheet part.
    parts <- entries[grepl("^xl/worksheets/sheet[0-9]+\\.xml$", entries)]
    if (length(parts) == 0) {
      stop("No worksheets found in ", basename(path), call. = FALSE)
    }
    sheet_names <- "Sheet1"
    targets <- parts[1]
  } else {
    targets <- vapply(sheet_rids, function(rid) {
      t <- rel_targ[match(rid, rel_ids)]
      if (is.na(t)) return(NA_character_)
      t <- sub("^/xl/", "", t)
      t <- sub("^xl/", "", t)
      if (t %in% entries) t else {
        hit <- entries[endsWith(tolower(entries), tolower(t))]
        if (length(hit)) hit[1] else NA_character_
      }
    }, character(1))
  }

  out <- list()
  for (i in seq_along(sheet_names)) {
    tgt <- targets[i]
    if (is.na(tgt)) next
    ws <- parse_worksheet(read_zip_entry(path, tgt), shared)
    out[[sheet_names[i]]] <- worksheet_to_df(ws)
  }
  if (length(out) == 0) {
    stop("No readable worksheet in ", basename(path), call. = FALSE)
  }
  out
}

# ---------------------------------------------------------------------------
#  Writing
# ---------------------------------------------------------------------------

.crc_table <- local({
  poly <- -306674912L          # 0xEDB88320 as a signed 32-bit integer
  t0 <- integer(256)
  for (i in 0:255) {
    cc <- i
    for (k in 1:8) {
      cc <- if (bitwAnd(cc, 1L) == 1L) bitwXor(poly, bitwShiftR(cc, 1L))
            else bitwShiftR(cc, 1L)
    }
    t0[i + 1L] <- cc
  }
  shift_xor <- function(t) bitwXor(bitwShiftR(t, 8L), t0[bitwAnd(t, 255L) + 1L])
  t1 <- shift_xor(t0); t2 <- shift_xor(t1); t3 <- shift_xor(t2)
  list(t0 = t0, t1 = t1, t2 = t2, t3 = t3)
})

#' CRC-32 (IEEE 802.3) of a raw vector
#'
#' Standard reflected CRC-32 (init 0xFFFFFFFF, final XOR 0xFFFFFFFF), using
#' slicing-by-4 so that four bytes are consumed per R-level iteration.
#' @keywords internal
crc32_raw <- function(bytes) {
  n <- length(bytes)
  if (n == 0L) return(0)
  v <- as.integer(bytes)
  T0 <- .crc_table$t0; T1 <- .crc_table$t1
  T2 <- .crc_table$t2; T3 <- .crc_table$t3

  crc <- -1L                                  # 0xFFFFFFFF as signed 32-bit
  i <- 1L
  lim <- n - 3L
  while (i <= lim) {
    b <- bitwOr(bitwOr(v[i], bitwShiftL(v[i + 1L], 8L)),
                bitwOr(bitwShiftL(v[i + 2L], 16L),
                       bitwShiftL(v[i + 3L], 24L)))
    x <- bitwXor(crc, b)
    crc <- bitwXor(
      bitwXor(T3[bitwAnd(x, 255L) + 1L],
              T2[bitwAnd(bitwShiftR(x, 8L), 255L) + 1L]),
      bitwXor(T1[bitwAnd(bitwShiftR(x, 16L), 255L) + 1L],
              T0[bitwAnd(bitwShiftR(x, 24L), 255L) + 1L]))
    i <- i + 4L
  }
  while (i <= n) {
    crc <- bitwXor(bitwShiftR(crc, 8L),
                   T0[bitwAnd(bitwXor(crc, v[i]), 255L) + 1L])
    i <- i + 1L
  }
  crc <- bitwXor(crc, -1L)
  if (crc < 0) crc + 4294967296 else crc
}

le16 <- function(v) as.raw(c(v %% 256L, (v %/% 256L) %% 256L))
le32 <- function(v) {
  v <- as.numeric(v)
  as.raw(c(v %% 256, floor(v / 256) %% 256,
           floor(v / 65536) %% 256, floor(v / 16777216) %% 256))
}

#' Write a named list of data.frames to an .xlsx file
#'
#' Dispatches to \pkg{writexl} when available (fast, but optional) and falls
#' back to the built-in writer otherwise. Force a specific engine with
#' \code{options(ADEPT.xlsx_engine = "internal")}.
#'
#' @param sheets Named list of data.frames.
#' @param path Output path.
#' @param engine \code{"auto"} (default), \code{"writexl"} or \code{"internal"}.
#' @keywords internal
write_xlsx_base <- function(sheets, path,
                            engine = getOption("ADEPT.xlsx_engine", "auto")) {
  stopifnot(length(sheets) > 0)
  engine <- match.arg(engine, c("auto", "writexl", "internal"))

  # writexl (zero R dependencies, pure C) is 10-100x faster when present.
  # It is only a suggestion: the built-in writer below is always available.
  use_writexl <- identical(engine, "writexl") ||
    (identical(engine, "auto") && requireNamespace("writexl", quietly = TRUE))

  if (use_writexl) {
    ok <- tryCatch({
      writexl::write_xlsx(sheets, path)
      TRUE
    }, error = function(e) FALSE)
    if (ok) return(invisible(path))
    warning("writexl failed; falling back to the built-in writer.",
            call. = FALSE)
  }
  write_xlsx_internal(sheets, path)
}

#' Write an .xlsx workbook using the built-in (dependency-free) writer
#'
#' Entries are stored uncompressed, so no deflate library is needed. The
#' resulting file is a valid OOXML workbook and is read correctly by Excel,
#' LibreOffice, readxl and \code{\link{adept_read}}.
#'
#' @param sheets Named list of data.frames.
#' @param path Output path.
#' @keywords internal
write_xlsx_internal <- function(sheets, path) {
  nm <- names(sheets)
  if (is.null(nm) || any(!nzchar(nm))) nm <- paste0("Sheet", seq_along(sheets))
  nm <- substr(gsub("[\\[\\]:*?/\\\\]", "_", nm), 1, 31)
  nm <- make.unique(nm)

  max_cols <- max(vapply(sheets, ncol, integer(1)), 0L)
  ref_letters <- if (max_cols > 0) .num_to_col_letters(seq_len(max_cols)) else character(0)

  cell_xml <- function(df, row_offset = 1L) {
    nr <- nrow(df); nc <- ncol(df)
    if (nc == 0L || nr == 0L) return("")
    cells <- vector("list", nc)
    abs_row <- seq_len(nr) + row_offset
    for (j in seq_len(nc)) {
      v <- df[[j]]
      letter <- ref_letters[j]          # pre-computed A, B, ..., AA
      if (is.logical(v)) v <- as.numeric(v)
      if (is.numeric(v)) {
        txt <- as.character(v)
        cells[[j]] <- ifelse(is.na(v),
                             paste0('<c r="', letter, abs_row, '"/>'),
                             paste0('<c r="', letter, abs_row, '"><v>',
                                    txt, "</v></c>"))
      } else {
        txt <- as.character(v)
        cells[[j]] <- ifelse(is.na(txt),
                             paste0('<c r="', letter, abs_row, '"/>'),
                             paste0('<c r="', letter, abs_row,
                                    '" t="inlineStr"><is><t>',
                                    xml_escape(txt), "</t></is></c>"))
      }
    }
    # Vectorised row assembly: paste the nc column-vectors element-wise.
    rowbody <- if (nc == 1L) cells[[1]] else do.call(paste0, cells)
    paste0('<row r="', abs_row, '">', rowbody, "</row>", collapse = "")
  }

  # ---- build parts -------------------------------------------------------
  sheet_parts <- character(length(sheets))
  for (i in seq_along(sheets)) {
    df <- sheets[[i]]
    hdr <- if (ncol(df) > 0) {
      paste0('<row r="1">',
             paste0('<c r="', ref_letters[seq_len(ncol(df))],
                    '1" t="inlineStr"><is><t>',
                    xml_escape(names(df)), "</t></is></c>", collapse = ""),
             "</row>")
    } else ""
    body <- if (nrow(df) > 0) cell_xml(df, row_offset = 1L) else ""
    sheet_parts[i] <- paste0(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
      "<sheetData>", hdr, body, "</sheetData></worksheet>"
    )
  }

  overrides <- paste0(vapply(seq_along(sheets), function(i)
    sprintf('<Override PartName="/xl/worksheets/sheet%d.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>',
            i), character(1)), collapse = "")

  parts <- list(
    "[Content_Types].xml" = paste0(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
      '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
      '<Default Extension="xml" ContentType="application/xml"/>',
      '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
      overrides, "</Types>"),

    "_rels/.rels" = paste0(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
      '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>',
      "</Relationships>"),

    "xl/workbook.xml" = paste0(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" ',
      'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>',
      paste0(sprintf('<sheet name="%s" sheetId="%d" r:id="rId%d"/>',
                     xml_escape(nm), seq_along(nm), seq_along(nm)),
             collapse = ""),
      "</sheets></workbook>"),

    "xl/_rels/workbook.xml.rels" = paste0(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
      paste0(sprintf('<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet%d.xml"/>',
                     seq_along(sheets), seq_along(sheets)), collapse = ""),
      "</Relationships>")
  )
  for (i in seq_along(sheets)) {
    parts[[sprintf("xl/worksheets/sheet%d.xml", i)]] <- sheet_parts[i]
  }

  .write_zip_stored(parts, path)
  invisible(path)
}

#' Column index to Excel letters
#' @keywords internal
.num_to_col_letters <- function(n) {
  out <- character(length(n))
  for (k in seq_along(n)) {
    v <- n[k]; s <- ""
    while (v > 0) {
      r <- (v - 1) %% 26
      s <- paste0(intToUtf8(65 + r), s)
      v <- (v - 1) %/% 26
    }
    out[k] <- s
  }
  out
}

#' Write a set of named text parts as a ZIP archive with stored entries
#'
#' @param parts Named character vector: archive path -> file content (UTF-8).
#' @param path Output .zip/.xlsx path.
#' @keywords internal
.write_zip_stored <- function(parts, path) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)

  central <- vector("list", length(parts))
  offset <- 0
  dos_time <- 0; dos_date <- 33      # 1980-01-01

  for (i in seq_along(parts)) {
    name <- names(parts)[i]
    nb   <- charToRaw(enc2utf8(name))
    data <- charToRaw(enc2utf8(parts[[i]]))
    crc  <- crc32_raw(data)
    sz   <- length(data)

    local <- c(le32(0x04034b50), le16(20L), le16(0L), le16(0L),
               le16(dos_time), le16(dos_date),
               le32(crc), le32(sz), le32(sz),
               le16(length(nb)), le16(0L), nb)
    writeBin(local, con)
    if (sz > 0) writeBin(data, con)

    central[[i]] <- c(le32(0x02014b50), le16(20L), le16(20L), le16(0L),
                      le16(0L), le16(dos_time), le16(dos_date),
                      le32(crc), le32(sz), le32(sz),
                      le16(length(nb)), le16(0L), le16(0L),
                      le16(0L), le16(0L), le32(0L), le32(offset), nb)
    offset <- offset + length(local) + sz
  }

  cd_offset <- offset
  cd_size <- 0
  for (i in seq_along(central)) {
    writeBin(central[[i]], con)
    cd_size <- cd_size + length(central[[i]])
  }
  writeBin(c(le32(0x06054b50), le16(0L), le16(0L),
             le16(length(parts)), le16(length(parts)),
             le32(cd_size), le32(cd_offset), le16(0L)), con)
  invisible(path)
}
