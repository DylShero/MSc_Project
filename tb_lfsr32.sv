// tb_lfsr32.sv  -- simulation 
// Drives lfsr32 and checks every output word against the baseline model.
//
// Verification steps:
//   1. gen_expected (from C++ baseline model) writes expected.hex
//   2. this testbench loads the same seed, advances the DUT, and compares
//      each word to expected.hex via $readmemh.


`timescale 1ns/1ps

module tb_lfsr32;
    localparam int N = 20;                 // words to check (== lines in expected.hex)

    logic        clk = 0, rst, en, load;
    logic [31:0] seed_in, rnd_out;
    logic [31:0] expected [0:N-1];
    integer      k, errors = 0;

    // Device under test
    lfsr32 dut (.clk(clk), .rst(rst), .en(en), .load(load),
                .seed_in(seed_in), .rnd_out(rnd_out));

    // 100 MHz clock
    always #5 clk = ~clk;

    initial begin
        $readmemh("expected.hex", expected); //Expected output answers in hexadeciaml

        // reset 
        rst = 1; en = 0; load = 0; seed_in = 32'h1234_5678;   // same seed as baseline
        @(negedge clk); @(negedge clk);

        // load the initial state
        rst = 0; load = 1;
        @(negedge clk);

        // advance and compare
        // The register holds the LFSR STATE. The seed we just loaded is a state, not an output word. 
        // the baseline model's first word is the state after one advance

        load = 0; en = 1;
        for (k = 0; k < N; k = k + 1) begin
            @(posedge clk); #1;            // let the non-blocking update settle
            if (rnd_out !== expected[k]) begin
                errors = errors + 1;
                $display("  MISMATCH k=%0d  hw=%08x  golden=%08x", k, rnd_out, expected[k]);
            end
        end

        if (errors == 0) $display(">>> PASS: all %0d words bit-exact with baseline <<<", N);
        else             $display(">>> FAIL: %0d mismatch(es) <<<", errors);
        $finish;
    end
endmodule
