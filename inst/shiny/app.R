# ===========================================================================
#  ADEPT 可视化工作台 (Shiny)
# ===========================================================================
#
#  运行方式（三选一）：
#    1) R 控制台：  ADEPT::adept_gui()
#    2) RStudio：   打开本文件 -> Run App
#    3) 命令行：    shiny::runApp(system.file("shiny", package = "ADEPT"))
#
#  想改界面直接改这个文件即可，不需要重新编译 / 重新安装 R 包。
#
#  可选依赖：
#    DT —— 交互式表格（未安装时自动降级为静态表格）
#
# ---------------------------------------------------------------------------

library(shiny)
library(ADEPT)

HAS_DT <- requireNamespace("DT", quietly = TRUE)
if (HAS_DT) library(DT)

# ---------------------------------------------------------------------------
#  小工具
# ---------------------------------------------------------------------------

# 统一的表格 UI：有 DT 用 DT，没有就退回基础 tableOutput
table_ui <- function(id, height = "460px") {
  if (HAS_DT) {
    DT::DTOutput(id, height = height)
  } else {
    div(style = sprintf("max-height:%s; overflow:auto;", height),
        tableOutput(id))
  }
}

# 统一的表格渲染。注意 output 必须显式传入 —— 这是 server 里的响应式对象，
# 在普通函数中不可见。
render_table <- function(output, id, df, digits = 4) {
  if (HAS_DT) {
    output[[id]] <- DT::renderDT({
      d <- df()
      req(d)
      num <- which(vapply(d, is.numeric, logical(1)))
      tbl <- DT::datatable(
        d,
        rownames = FALSE,
        extensions = "Buttons",
        options = list(pageLength = 25, scrollX = TRUE,
                       dom = "Bfrtip",
                       buttons = c("copy", "csv", "excel"))
      )
      if (length(num) > 0) {
        tbl <- DT::formatRound(tbl, columns = num, digits = digits)
      }
      tbl
    }, server = TRUE)
  } else {
    output[[id]] <- renderTable({ req(df()); df() },
                                digits = digits, spacing = "s")
  }
}

# 把 data.frame 渲染成 Bootstrap 表格（不引入 markdown 等额外依赖）
md_table <- function(df) {
  tags$table(
    class = "table table-condensed table-striped",
    tags$thead(tags$tr(lapply(colnames(df), tags$th))),
    tags$tbody(lapply(seq_len(nrow(df)), function(i) {
      tags$tr(lapply(df[i, ], tags$td))
    }))
  )
}

help_ui <- function() {
  tagList(
    tags$h4("1. 输入数据"),
    tags$p("支持三种 Excel 格式，列名必须完全一致（注意大小写）："),
    md_table(data.frame(
      `格式` = c("1 · 直接年龄", "2 · 原始计数", "3 · 同位素比值"),
      `需要的列` = c("Analysis, Time, Age68, Age75, Age76",
                     "Analysis, Time, Pb206, Pb207, U238",
                     "Analysis, Time, Pb206_U238, Pb207_U235, Pb207_Pb206"),
      check.names = FALSE)),
    tags$p("格式 2 / 3 需要额外安装 ", tags$code("IsoplotR"),
           " 来把比值换算成年龄：",
           tags$code("install.packages(\"IsoplotR\")")),
    tags$p("列名里的空格与斜杠会自动替换为下划线。除上表之外的",
           tags$strong("数值列"),"（如微量元素）会被自动识别，",
           "并在结果中输出其在各年龄坪内的均值。"),

    tags$h4("2. 参数含义"),
    md_table(data.frame(
      `参数` = c("剥蚀时间窗口", "方差阈值", "最小坪宽 (s)",
                 "年龄上/下限", "筛选方向", "MCMC"),
      `说明` = c("只保留这段时间内的信号，避开开头爬升与结尾穿透污染（默认 29–58 s）",
                 "年龄坪内部归一化方差上限。调小 → 更严格、坪更少；调大 → 更宽松（默认 0.1192）",
                 "短于该时长的坪被剔除（第一段豁免，因为信号起点常被截断）",
                 "落在区间外的坪直接剔除",
                 "Forward = 年龄随深度递增；Reverse = 年龄随深度递减",
                 "贝叶斯变点后验估计，需 mcp + JAGS，未安装时该列留空"),
      check.names = FALSE)),

    tags$h4("3. 建议工作流"),
    tags$ol(
      tags$li("用默认参数 + 「载入自带示例数据」跑一遍，确认流程通。"),
      tags$li("在「数据预览」页确认输入格式被正确识别。"),
      tags$li("在「深度剖面」页逐个锆石看图，判断年龄趋势是递增还是递减，据此选择方向。"),
      tags$li("坪太少 → 调大方差阈值或调小最小坪宽；坪太碎 → 反向调整。"),
      tags$li("参数定下来后，在论文方法部分写明所用参数。")
    ),

    tags$h4("4. 注意"),
    tags$p(tags$strong("年龄坪是算法在给定参数下的最优分段，"),
           "并不等同于地质上真实的结晶事件。对铀矿、热液蚀变锆石等体系，",
           "U-Pb 体系可能遭受多期扰动，此时年龄坪未必具有地质意义，",
           "务必结合 CL / BSE 图像与微量元素剖面综合判断。")
  )
}
# 画一个深度剖面（年龄坪用红色横线 + 不确定性带表示）
# 用 base graphics：包本体已不再依赖 ggplot2，界面也不需要它。
draw_profile <- function(profile, cex_point = 0.6) {
  df  <- profile$data
  seg <- profile$segments

  op <- graphics::par(mar = c(4.0, 4.3, 2.4, 1.2),
                      mgp = c(2.5, 0.6, 0), tcl = -0.3,
                      xaxs = "r", yaxs = "r",
                      cex.axis = 0.95, cex.lab = 1.05, cex.main = 1.1)
  on.exit(graphics::par(op), add = TRUE)

  graphics::plot(df$Time, df$Raw_Age, type = "n",
                 xlab = "Ablation time (s)", ylab = "Age (Ma)",
                 main = paste("Analysis:", profile$analysis))

  keep <- !is.na(seg$Filter_4)
  if (any(keep)) {
    graphics::rect(
      seg$Start[keep],
      seg$Segment_Mean[keep] - seg$Total_uncertainty[keep],
      seg$End[keep],
      seg$Segment_Mean[keep] + seg$Total_uncertainty[keep],
      col = grDevices::adjustcolor("red", 0.15), border = NA)
  }
  graphics::points(df$Time, df$Raw_Age, pch = 16, cex = cex_point,
                   col = grDevices::adjustcolor("#f46f20", 0.4))
  graphics::lines(df$Time, df$loess_Age, lwd = 1.6)
  if (any(keep)) {
    graphics::segments(seg$Start[keep], seg$Segment_Mean[keep],
                       seg$End[keep],   seg$Segment_Mean[keep],
                       col = "red", lwd = 2)
  }
  graphics::box()
  invisible(profile)
}

# ---------------------------------------------------------------------------
#  UI
# ---------------------------------------------------------------------------

ui <- fluidPage(

  titlePanel("ADEPT · 锆石 U-Pb 深度剖面年龄坪自动提取"),

  sidebarLayout(

    sidebarPanel(
      width = 3,

      fileInput("file", "选择输入 Excel (.xlsx)",
                accept = c(".xlsx", ".xls")),
      actionButton("use_demo", "载入自带示例数据",
                   class = "btn-default btn-sm", width = "100%"),

      tags$hr(),
      tags$h5("剥蚀时间窗口"),
      numericInput("lower_time", "起始 (s)", value = 29, min = 0, step = 1),
      numericInput("upper_time", "结束 (s)", value = 58, min = 0, step = 1),

      tags$hr(),
      tags$h5("年龄坪判据"),
      numericInput("var_thr", "方差阈值", value = 0.1192, min = 0, step = 0.001),
      numericInput("min_res", "最小坪宽 (s)", value = 5, min = 0, step = 1),
      numericInput("min_age", "年龄下限 (Ma)", value = 0, step = 1),
      numericInput("max_age", "年龄上限 (Ma)", value = 4540, step = 1),
      selectInput("direction", "筛选方向",
                  choices = c("Forward（年龄递增）" = "Forward",
                              "Reverse（年龄递减）" = "Reverse"),
                  selected = "Forward"),

      tags$hr(),
      tags$h5("其它"),
      numericInput("chunk", "每个锆石的行数", value = 411, min = 1, step = 1),
      checkboxInput("mcmc", "运行贝叶斯 MCMC（需 mcp + JAGS）", value = FALSE),
      checkboxInput("make_plots", "生成深度剖面图", value = TRUE),

      tags$hr(),
      actionButton("run", "开始处理", class = "btn-primary", width = "100%"),
      tags$br(), tags$br(),
      downloadButton("dl_xlsx", "下载 Excel 结果", width = "100%"),
      tags$br(),
      downloadButton("dl_pdf", "下载全部剖面图 (PDF)", width = "100%")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("说明", help_ui()),

        tabPanel("数据预览",
                 tags$h4("输入文件结构"),
                 verbatimTextOutput("file_info"),
                 tags$h4("前 12 行"),
                 table_ui("tbl_preview", height = "340px")),

        tabPanel("深度剖面",
                 fluidRow(
                   column(6, selectInput("pick_profile", "选择锆石",
                                         choices = NULL)),
                   column(3, numericInput("pw", "图宽 (cm)",
                                          value = 14, min = 6, step = 1)),
                   column(3, numericInput("ph", "图高 (cm)",
                                          value = 8, min = 4, step = 1))
                 ),
                 plotOutput("profile_plot", height = "480px"),
                 downloadButton("dl_one", "下载当前图 (PDF)")),

        tabPanel("年龄坪汇总",
                 tags$p(tags$b("仅包含通过全部四步筛选的年龄坪。")),
                 table_ui("tbl_summary")),

        tabPanel("完整结果",
                 tags$p("全部年龄坪（含被筛掉的）与全部统计量。"),
                 table_ui("tbl_full")),

        tabPanel("运行日志", verbatimTextOutput("log_out"))
      )
    )
  )
)

# ---------------------------------------------------------------------------
#  Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  rv <- reactiveValues(result   = NULL,
                       profiles = NULL,
                       xlsx     = NULL,
                       logs     = character(0))

  src_path <- reactiveVal(NULL)

  reset_result <- function() {
    rv$result <- NULL
    rv$profiles <- NULL
    rv$xlsx <- NULL
  }

  observeEvent(input$file, {
    req(input$file)
    src_path(input$file$datapath)
    reset_result()
    rv$logs <- paste0("已载入文件: ", input$file$name)
  })

  observeEvent(input$use_demo, {
    demo <- system.file("extdata", "Input(1sample).xlsx", package = "ADEPT")
    if (nzchar(demo)) {
      src_path(demo)
      reset_result()
      rv$logs <- paste0("已载入示例数据: ", basename(demo))
    } else {
      showNotification("未找到自带示例数据。", type = "warning")
    }
  })

  # ---- 数据预览 ----------------------------------------------------------
  # 用包自带的 xlsx 读取器，不需要 readxl / openxlsx
  preview <- reactive({
    p <- src_path()
    req(p)
    sheets <- try(ADEPT::adept_read(p), silent = TRUE)
    if (inherits(sheets, "try-error") || length(sheets) == 0) return(NULL)
    d <- as.data.frame(sheets[[1]])
    list(sheets = names(sheets), colnames = colnames(d),
         head = utils::head(d, 12))
  })

  output$file_info <- renderText({
    pv <- preview()
    req(pv)
    cn <- pv$colnames
    fmt <- if (all(c("Age68", "Age75", "Age76") %in% cn)) {
      "格式 1 —— 直接年龄 (Age68 / Age75 / Age76)"
    } else if (all(c("Pb206", "Pb207", "U238") %in% cn)) {
      "格式 2 —— 原始同位素计数（需 IsoplotR 换算年龄）"
    } else if (all(c("Pb206_U238", "Pb207_U235", "Pb207_Pb206") %in% cn)) {
      "格式 3 —— 同位素比值（需 IsoplotR 换算年龄）"
    } else {
      "无法识别输入格式，请检查列名"
    }
    paste0("工作表 (", length(pv$sheets), "): ",
           paste(pv$sheets, collapse = ", "), "\n",
           "列数: ", length(cn), "\n",
           "列名: ", paste(cn, collapse = ", "), "\n\n",
           "识别结果: ", fmt)
  })

  render_table(output, "tbl_preview", reactive({
    pv <- preview(); req(pv); pv$head
  }))

  # ---- 运行 --------------------------------------------------------------
  observeEvent(input$run, {
    p <- src_path()
    if (is.null(p)) {
      showNotification("请先选择输入文件，或载入自带示例数据。",
                       type = "warning")
      return(NULL)
    }

    out_xlsx <- tempfile(fileext = "_ADEPT_Output.xlsx")
    logs <- character(0)

    res <- try(withProgress(message = "正在处理，请稍候…", value = 0, {
      ADEPT::adept(
        file_path              = p,
        chunk_size             = input$chunk,
        lower_ablation_time    = input$lower_time,
        upper_ablation_time    = input$upper_time,
        max_age_limit          = input$max_age,
        min_age_limit          = input$min_age,
        min_plateau_resolution = input$min_res,
        variance_threshold     = input$var_thr,
        filter_direction       = input$direction,
        mcmc                   = input$mcmc,
        make_plots             = input$make_plots,
        save_plots_to_disk     = FALSE,
        output_path            = out_xlsx,
        keep_profiles          = TRUE,
        verbose                = FALSE,
        progress = function(fraction, detail) {
          logs <<- c(logs, sprintf("[%3.0f%%] %s", fraction * 100, detail))
          incProgress(amount = 0, detail = detail)
        }
      )
    }), silent = TRUE)

    if (inherits(res, "try-error")) {
      msg <- conditionMessage(attr(res, "condition"))
      rv$logs <- c(logs, paste("ERROR:", msg))
      showNotification(paste("处理失败：", msg), type = "error", duration = 10)
      return(NULL)
    }

    rv$logs     <- logs
    rv$xlsx     <- out_xlsx
    rv$result   <- res
    rv$profiles <- res$profiles

    if (!is.null(res$profiles) && length(res$profiles) > 0) {
      labels <- vapply(seq_along(res$profiles), function(i) {
        pr <- res$profiles[[i]]
        sprintf("%02d | %s | %s", i, pr$analysis, pr$sheet)
      }, character(1))
      updateSelectInput(session, "pick_profile",
                        choices = setNames(as.list(seq_along(res$profiles)),
                                           labels),
                        selected = 1)
    }

    n_plateau <- if (!is.null(res$summary)) nrow(res$summary) else 0
    showNotification(
      sprintf("完成：%d 个锆石，识别到 %d 个年龄坪",
              length(res$profiles), n_plateau),
      type = "message")
    updateTabsetPanel(session, "tabs", selected = "深度剖面")
  })

  # ---- 剖面图 ------------------------------------------------------------
  current_profile <- reactive({
    req(rv$profiles)
    req(input$pick_profile)
    i <- suppressWarnings(as.integer(input$pick_profile))
    if (is.na(i) || i < 1 || i > length(rv$profiles)) return(NULL)
    rv$profiles[[i]]
  })

  output$profile_plot <- renderPlot({
    pr <- current_profile()
    req(pr)
    draw_profile(pr)
  })

  output$dl_one <- downloadHandler(
    filename = function() {
      pr <- current_profile()
      nm <- if (is.null(pr)) "profile" else make.names(pr$analysis)
      paste0("ADEPT_", nm, ".pdf")
    },
    content = function(file) {
      pr <- current_profile()
      req(pr)
      grDevices::pdf(file, width = input$pw / 2.54, height = input$ph / 2.54)
      draw_profile(pr)
      grDevices::dev.off()
    },
    contentType = "application/pdf"
  )

  output$dl_pdf <- downloadHandler(
    filename = "ADEPT_all_profiles.pdf",
    content = function(file) {
      req(rv$profiles)
      grDevices::pdf(file, width = input$pw / 2.54, height = input$ph / 2.54)
      for (pr in rv$profiles) draw_profile(pr)
      grDevices::dev.off()
    },
    contentType = "application/pdf"
  )

  output$dl_xlsx <- downloadHandler(
    filename = function() "ADEPT_results.xlsx",
    content = function(file) {
      req(rv$xlsx)
      file.copy(rv$xlsx, file, overwrite = TRUE)
    },
    contentType =
      "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )

  # ---- 结果表 ------------------------------------------------------------
  render_table(output, "tbl_summary", reactive({
    req(rv$result); rv$result$summary
  }))

  render_table(output, "tbl_full", reactive({
    req(rv$result); rv$result$full
  }))

  output$log_out <- renderText({
    if (length(rv$logs) == 0) return("（尚无日志）")
    paste(rv$logs, collapse = "\n")
  })
}

shinyApp(ui, server)
