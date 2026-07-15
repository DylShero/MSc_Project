`timescale 1ns/1ps
module tb_xorshift128;
    localparam int N = 20;
    logic         clk = 0, rst, en, load;
    logic [127:0] seed_in;
    logic [31:0]  rnd_out;
    logic [31:0]  expected [0:N-1];
    integer       k, errors = 0;

    xorshift128 dut (.clk(clk), .rst(rst), .en(en), .load(load),
                     .seed_in(seed_in), .rnd_out(rnd_out));

    always #5 clk = ~clk; //Delay by 5ns

    initial begin
        $readmemh("expected_xorshift.hex", expected);
        rst = 1; en = 0; load = 0;
        seed_in = 128'h12345678_9ABCDEF0_0FEDCBA9_87654321;  // {x,y,z,w}, same as baseline
        @(negedge clk); @(negedge clk);
        rst = 0; load = 1; @(negedge clk);
        load = 0; en = 1;
        for (k = 0; k < N; k = k + 1) begin
            @(posedge clk); #1;
            if (rnd_out !== expected[k]) begin
                errors = errors + 1;
                $display("  MISMATCH k=%0d  hw=%08x  baseline=%08x", k, rnd_out, expected[k]);
            end
        end
        if (errors == 0) $display(">>> PASS: all %0d words match baseline <<<", N);
        else             $display(">>> FAIL: %0d mismatch(es) <<<", errors);
        $finish;
    end
endmodule
