test_that("result lists retain zero, one and multiple items across R and Python", {
  py <- batchdx_python()
  if (is.null(py)) skip("no Python interpreter can import batchdx.cyto")
  for (values in list(character(0), "node 12", c("node 12", "node 34"))) {
    path <- file.path(batchdx_workdir("result"), "result.json")
    batchdx_write_result(path, list(status = "ok", warnings = values,
                                   confounded_nodes = values))
    script <- paste(
      "import json, sys",
      "with open(sys.argv[1]) as f: payload = json.load(f)",
      "assert payload['status'] == 'ok'",
      "assert isinstance(payload['warnings'], list)",
      "assert isinstance(payload['confounded_nodes'], list)",
      "print(json.dumps(payload))", sep = "\n"
    )
    out <- system2(py, c("-c", shQuote(script), shQuote(path)), stdout = TRUE, stderr = TRUE)
    expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
    payload <- jsonlite::fromJSON(paste(out, collapse = "\n"), simplifyVector = FALSE)
    expect_identical(payload$warnings, as.list(values))
    expect_identical(payload$confounded_nodes, as.list(values))
  }
})
