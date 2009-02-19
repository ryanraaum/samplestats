require 'rake/clean'

#
# Platform specific options
#
if RUBY_PLATFORM =~ /win32/i then
  COMPILE             = "cl"
  COMPILE_FLAGS       = "/O2 /W3 /c /D_CRT_SECURE_NO_DEPRECATE /nologo"
  COMPILE_OBJECT_FLAG = "/Fo"
  LINK                = "link"
  LINK_FLAGS          = "/nologo /out:"
  
  EXEC_EXTENSION      = ".exe"
  EXEC_PREFIX         = ""
else
  COMPILE             = "gcc"
  COMPILE_FLAGS       = "-O2 -c"
  COMPILE_OBJECT_FLAG = "-o "
  LINK                = "gcc"
  LINK_FLAGS          = "-o "
  
  EXEC_EXTENSION      = ""
  EXEC_PREFIX         = "./"
end

#
# Executable names
#
TESTGETOPTPROG        = 'test_simple_getopt'  + EXEC_EXTENSION
TESTUNICFREQSPROG     = 'test_unic_freqs'     + EXEC_EXTENSION
SAMPLESTATSPROG       = 'sample_stats'        + EXEC_EXTENSION
SAMPLESTATSPROG2      = 'sample_stats2'       + EXEC_EXTENSION
SAMPLESTATSPROG3      = 'sample_stats3'       + EXEC_EXTENSION

#
# Some lists to be used in clean and clobber tasks
#
EXECUTABLES           = [ TESTGETOPTPROG, 
                          SAMPLESTATSPROG, 
                          SAMPLESTATSPROG2,
                          SAMPLESTATSPROG3 ]
SRC                   = FileList['*.c']
OBJ                   = SRC.collect { |fn| File.basename(fn).ext('o') }
TESTFILES             = [ "small_theta_ms_output",
                          "small_theta_sample_stats_output",
                          "small_theta_sample_stats2_output",
                          "big_theta_ms_output",
                          "big_theta_sample_stats_output",
                          "big_theta_sample_stats2_output",
                          "ss2_out",
                          "ss3_out",
                          "opttestout",
                          "treefile",
                          "onesmallseqgen", 
                          "manysmallseqgen",
                          "onebigseqgen", 
                          "manybigseqgen" ]
TEMPFILES             = [ "seedms" ]

#
# Configure rake-supplied clean and clobber tasks
#
CLEAN.include(OBJ, TESTFILES, TEMPFILES)
CLOBBER.include(EXECUTABLES)

# Default rule for object files without dependencies
rule '.o' => ['.c'] do |t|
    sh "#{COMPILE} #{COMPILE_FLAGS} #{t.prerequisites.join(' ')} #{COMPILE_OBJECT_FLAG}#{t.name}"
end

#
# Rules to build executables
#

file TESTGETOPTPROG => ["test_simple_getopt.o", "simple_getopt.o"] do |t|
  sh "#{LINK} #{LINK_FLAGS}#{t.name} #{t.prerequisites.join(' ')}"
end

file TESTUNICFREQSPROG => ["test_unic_freqs.o", "r2.o" ] do |t|
  sh "#{LINK} #{LINK_FLAGS}#{t.name} #{t.prerequisites.join(' ')}"
end

file SAMPLESTATSPROG => ["sample_stats.o", "tajd.o"] do |t|
  sh "#{LINK} #{LINK_FLAGS}#{t.name} #{t.prerequisites.join(' ')}"
end

file SAMPLESTATSPROG2 => ["sample_stats2.o", "tajd.o", "fs.o", "r2.o", "simple_getopt.o"] do |t|
  sh "#{LINK} #{LINK_FLAGS}#{t.name} #{t.prerequisites.join(' ')}"
end

file SAMPLESTATSPROG3 => ["sample_stats3.o", "tajd.o", "r2.o", "simple_getopt.o"] do |t|
  sh "#{LINK} #{LINK_FLAGS}#{t.name} #{t.prerequisites.join(' ')}"
end

#
# Tasks to build excutables (and groups of executables)
#

desc "build Hudson's sample_stats"
task :build_sample_stats  => [SAMPLESTATSPROG]

desc "build Raaums's sample_stats2 for ms output"
task :build_sample_stats2 => [SAMPLESTATSPROG2]

desc "build Raaums's sample_stats3 for seq-gen output"
task :build_sample_stats3 => [SAMPLESTATSPROG3]

task :test_unic_freqs => [TESTUNICFREQSPROG] do
  sh "#{TESTUNICFREQSPROG}"
end

#
# Rules to build some files used in tests
# These assume that Hudson's `ms` and Someone's `seq-gen` are in the PATH
#
file "small_theta_ms_output" do |t|
  sh "ms 100 1 -t 1 > #{t.name}" #, :verbose => false
end

file "big_theta_ms_output" do |t|
  sh "ms 100 1 -t 200 > #{t.name}" #, :verbose => false
end

file "onesmallseqgen" do |t|
  sh "ms 10 1 -T > treefile"
  sh "seq-gen -mHKY -l 40 < treefile > #{t.name}"
end
  
file "manysmallseqgen" do |t|
  sh "ms 10 5 -T > treefile"
  sh "seq-gen -mHKY -l 40 < treefile > #{t.name}"
end

file "onebigseqgen" do |t|
  sh "ms 10 1 -T > treefile"
  sh "seq-gen -mHKY -l 2000 < treefile > #{t.name}"
end
  
file "manybigseqgen" do |t|
  sh "ms 10 5 -T > treefile"
  sh "seq-gen -mHKY -l 2000 < treefile > #{t.name}"
end

#
# Tests
#
namespace :test do

  #
  # Some assertions that are not otherwise available, we're not
  # testing ruby code here...
  #
  def assert_passes(*msg)
    raise "Assertion failed! #{msg}" unless yield
  end
    
  def assert_fails(*msg)
    begin
      yield
    rescue
      return true
    end
    raise "Assertion failed to fail! #{msg}"
  end
  
  def assert_equal( a, b, msg="")
    raise "Assertion failed! #{a} is not equal to #{b}, #{msg}" unless a == b
  end
  
  #
  # Make sure the original sample stats runs without crashing or otherwise throwing an error
  #
  desc "test sample_stats"
  task :ss => [SAMPLESTATSPROG, "small_theta_ms_output", "big_theta_ms_output"] do
    puts ""
    puts "TEST: Running some ms output through sample_stats."
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG} < small_theta_ms_output > small_theta_sample_stats_output", :verbose => false
    end
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG} < big_theta_ms_output > big_theta_sample_stats_output", :verbose => false
    end
    puts "SUCCESS. Ran without problems."
  end
  
  #
  # Make sure that sample_stats2 produces output identical to the original sample_stats
  # in the default invocation on the same dataset
  #
  desc "test sample_stats2 vs. sample_stats"
  task :ss2vss => [SAMPLESTATSPROG, SAMPLESTATSPROG2, "small_theta_ms_output", "big_theta_ms_output"] do
    puts ""
    puts "TEST: Running some ms output through sample_stats2,\n      comparing it to sample_stats output."
    puts ""
    
    # generate statistics using both sample_stats and sample_stats2
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG} < small_theta_ms_output > small_theta_sample_stats_output", :verbose => false
    end
    assert_passes do 
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG} < big_theta_ms_output > big_theta_sample_stats_output", :verbose => false
    end
    assert_passes do 
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} < small_theta_ms_output > small_theta_sample_stats2_output", :verbose => false
    end
    assert_passes do 
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} < big_theta_ms_output > big_theta_sample_stats2_output", :verbose => false
    end
    
    # read in the results from the small theta set and compare
    ss  = (File.open("small_theta_sample_stats_output", "r") { |f| f.gets }).split(' ')
    ss2 = (File.open("small_theta_sample_stats2_output", "r") { |f| f.gets }).split(' ')
    (0...ss.length).each { |i| assert_equal(ss[i], ss2[i]) }
    
    # read in the results from the big theta set and compare
    ss  = (File.open("big_theta_sample_stats_output", "r") { |f| f.gets }).split(' ')
    ss2 = (File.open("big_theta_sample_stats2_output", "r") { |f| f.gets }).split(' ')
    (0...ss.length).each { |i| assert_equal(ss[i], ss2[i]) }
    
    puts "SUCCESS. Everything matches."
  end
  
  #
  # Make sure that sample_stats2 flags all work and produce output whose
  # format matches what is expected
  #
  desc "test sample_stats2 flags"
  task :ss2f => [SAMPLESTATSPROG2, "big_theta_ms_output"] do
    puts ""
    puts "Running tests of sample_stats2 flags."
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -S < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "ss:", val[0] )
    assert_passes { val[1] =~ /[0-9]+/ }
    assert_passes { val[1].to_i >= 0 }
    puts "-S (Segregating Sites)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -p < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "pi:", val[0] )
    assert_passes { val[1] =~ /[0-9]+\.[0-9]+/ }
    assert_passes { val[1].to_f >= 0.0 }
    puts "-p (pi, nucleotide diversity)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -F < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "thetaH:", val[0] )
    assert_passes { val[1] =~ /[0-9]+\.[0-9]+/ }
    assert_passes { val[1].to_f >= 0.0 }
    puts "-F (Fay's H)".ljust(40) + "OK"

    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -d < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "H:", val[0] )
    assert_passes { val[1] =~ /-?[0-9]+\.[0-9]+/ }
    puts "-d (H = pi - Fay's H)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -W < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "thetaW:", val[0] )
    assert_passes { val[1] =~ /[0-9]+\.[0-9]+/ }
    assert_passes { val[1].to_f >= 0.0 }
    puts "-W (Watterson's theta)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -D < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "D:", val[0] )
    assert_passes { val[1] =~ /-?[0-9]+\.[0-9]+/ }
    puts "-D (Tajima's D)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -H < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "homozygosity:", val[0] )
    assert_passes { val[1] =~ /[0-9]+\.[0-9]+/ }
    assert_passes { val[1].to_f <= 1.0 and val[1].to_f >= 0.0 }
    puts "-H (homozygosity)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -n < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "num_haplotypes:", val[0] )
    assert_passes { val[1] =~ /[0-9]+/ }
    assert_passes { val[1].to_i >= 0 }
    puts "-n (number of haplotypes)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -s < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "num_singletons:", val[0] )
    assert_passes { val[1] =~ /[0-9]+/ }
    assert_passes { val[1].to_i >= 0 }
    puts "-s (number of singletons)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -N < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "nss:", val[0] )
    assert_passes { val[1] =~ /[0-9]+/ }
    assert_passes { val[1].to_i >= 0 }
    puts "-N (number of singleton sites)".ljust(40) + "OK"

    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -f < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "hf:", val[0] )
    assert_passes { val[1] =~ /-?[0-9]+\.[0-9]+/ }
    puts "-f (mean haplotype frequency)".ljust(40) + "OK"

    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -i < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "ih:", val[0] )
    assert_passes { val[1] =~ /[0-9]+/ }
    assert_passes { val[1].to_i >= 0 }
    puts "-i (max number of identical haplotypes)".ljust(40) + "OK"

    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -R < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "r2:", val[0] )
    assert_passes { val[1] =~ /-?[0-9]+\.[0-9]+/ }
    puts "-R (Ramos-Onsins and Rozas' R2)".ljust(40) + "OK"

    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG2} -U < big_theta_ms_output > ss2_out", :verbose => false
    end
    val = (File.open("ss2_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "Fs:", val[0] )
    assert_passes { val[1] =~ /-?[0-9]+\.[0-9]+/ }
    puts "-U (Fu's Fs)".ljust(40) + "OK"

    puts "SUCCESS."
  end
  
  #
  # Make sure that sample_stats3 can read data
  #
  desc "test sample_stats3 ability to read data"
  task :ss3read => [SAMPLESTATSPROG3, "onesmallseqgen", "manysmallseqgen", "onebigseqgen", "manybigseqgen"] do 
    puts ""
    puts "Running tests of sample_stats3 data reading."
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} < onesmallseqgen"
    end
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} < manysmallseqgen"
    end
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} < onebigseqgen"
    end
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} < manybigseqgen"
    end
    
  end
  
  #
  # Make sure that sample_stats3 flags all work and produce output whose
  # format matches what is expected
  #
  desc "test sample_stats2 flags"
  task :ss3f => [SAMPLESTATSPROG3, "onebigseqgen"] do
    puts ""
    puts "Running tests of sample_stats3 flags."
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -S < onebigseqgen > ss3_out", :verbose => false
    end
    val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "ss:", val[0] )
    assert_passes { val[1] =~ /[0-9]+/ }
    assert_passes { val[1].to_i >= 0 }
    puts "-S (Segregating Sites)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -p < onebigseqgen > ss3_out", :verbose => false
    end
    val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "pi:", val[0] )
    assert_passes { val[1] =~ /[0-9]+\.[0-9]+/ }
    assert_passes { val[1].to_f >= 0.0 }
    puts "-p (pi, nucleotide diversity)".ljust(40) + "OK"
    
   # assert_passes do
   #   sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -F < onebigseqgen > ss3_out", :verbose => false
   # end
   # val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
   # assert_equal( "thetaH:", val[0] )
   # assert_passes { val[1] =~ /[0-9]+\.[0-9]+/ }
   # assert_passes { val[1].to_f >= 0.0 }
   # puts "-F (Fay's H)".ljust(40) + "OK"

   # assert_passes do
   #   sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -d < onebigseqgen > ss3_out", :verbose => false
   # end
   # val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
   # assert_equal( "H:", val[0] )
   # assert_passes { val[1] =~ /-?[0-9]+\.[0-9]+/ }
   # puts "-d (H = pi - Fay's H)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -W < onebigseqgen > ss3_out", :verbose => false
    end
    val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "thetaW:", val[0] )
    assert_passes { val[1] =~ /[0-9]+\.[0-9]+/ }
    assert_passes { val[1].to_f >= 0.0 }
    puts "-W (Watterson's theta)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -D < onebigseqgen > ss3_out", :verbose => false
    end
    val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "D:", val[0] )
    assert_passes { val[1] =~ /-?[0-9]+\.[0-9]+/ }
    puts "-D (Tajima's D)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -H < onebigseqgen > ss3_out", :verbose => false
    end
    val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "homozygosity:", val[0] )
    assert_passes { val[1] =~ /[0-9]+\.[0-9]+/ }
    assert_passes { val[1].to_f <= 1.0 and val[1].to_f >= 0.0 }
    puts "-H (homozygosity)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -n < onebigseqgen > ss3_out", :verbose => false
    end
    val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "num_haplotypes:", val[0] )
    assert_passes { val[1] =~ /[0-9]+/ }
    assert_passes { val[1].to_i >= 0 }
    puts "-n (number of haplotypes)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -s < onebigseqgen > ss3_out", :verbose => false
    end
    val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "num_singletons:", val[0] )
    assert_passes { val[1] =~ /[0-9]+/ }
    assert_passes { val[1].to_i >= 0 }
    puts "-s (number of singletons)".ljust(40) + "OK"
    
    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -N < onebigseqgen > ss3_out", :verbose => false
    end
    val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "nss:", val[0] )
    assert_passes { val[1] =~ /[0-9]+/ }
    assert_passes { val[1].to_i >= 0 }
    puts "-N (number of singleton sites)".ljust(40) + "OK"

   # assert_passes do
   #   sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -f < onebigseqgen > ss3_out", :verbose => false
   # end
   # val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
   # assert_equal( "hf:", val[0] )
   # assert_passes { val[1] =~ /-?[0-9]+\.[0-9]+/ }
   # puts "-f (mean haplotype frequency)".ljust(40) + "OK"

   # assert_passes do
   #   sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -i < onebigseqgen > ss3_out", :verbose => false
   # end
   # val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
   # assert_equal( "ih:", val[0] )
   # assert_passes { val[1] =~ /[0-9]+/ }
   # assert_passes { val[1].to_i >= 0 }
   # puts "-i (max number of identical haplotypes)".ljust(40) + "OK"

    assert_passes do
      sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -R < onebigseqgen > ss3_out", :verbose => false
    end
    val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
    assert_equal( "r2:", val[0] )
    assert_passes { val[1] =~ /-?[0-9]+\.[0-9]+/ }
    puts "-R (Ramos-Onsins and Rozas' R2)".ljust(40) + "OK"

   # assert_passes do
   #   sh "#{EXEC_PREFIX}#{SAMPLESTATSPROG3} -U < onebigseqgen > ss3_out", :verbose => false
   # end
   # val = (File.open("ss3_out", "r") { |infile| infile.gets }).split(' ')
   # assert_equal( "Fs:", val[0] )
   # assert_passes { val[1] =~ /-?[0-9]+\.[0-9]+/ }
   # puts "-U (Fu's Fs)".ljust(40) + "OK"

    puts "SUCCESS."
  end
  
  #
  # Make sure that simple_getopt works, catching options properly and throwing errors when expected
  #
  desc "test simple getopt"
  task :getopt => [TESTGETOPTPROG] do
    puts ""
    puts "Running tests of simple_getopt."
    puts ""

    puts "TEST: running '#{TESTGETOPTPROG}' with no parameters"
    assert_passes { sh("#{EXEC_PREFIX}#{TESTGETOPTPROG}", :verbose => false) }
    puts "SUCCESS."
    puts ""

    puts "TEST: running '#{TESTGETOPTPROG}' with known single no-parameter option"
    assert_passes { sh("#{EXEC_PREFIX}#{TESTGETOPTPROG} -b", :verbose => false) }
    puts "SUCCESS."
    puts ""
    
    puts "TEST: running '#{TESTGETOPTPROG}' with known single parameter option\n  (no space)"
    assert_passes { sh("#{EXEC_PREFIX}#{TESTGETOPTPROG} -f4", :verbose => false) }
    puts "SUCCESS."
    puts ""
    
    puts "TEST: running '#{TESTGETOPTPROG}' with known single parameter option\n  (with space)"
    assert_passes { sh("#{EXEC_PREFIX}#{TESTGETOPTPROG} -f 4", :verbose => false) }
    puts "SUCCESS."
    puts ""
    
    puts "TEST: running '#{TESTGETOPTPROG}' with known single parameter option\n  (without parameter!)"
    assert_fails { sh("#{EXEC_PREFIX}#{TESTGETOPTPROG} -f", :verbose => false) }
    puts "SUCCESS. Failed when it was supposed to."
    puts ""
    
    puts "TEST: running '#{TESTGETOPTPROG}' with long-style options\n  (not supported here)"
    assert_fails { sh("#{EXEC_PREFIX}#{TESTGETOPTPROG} --version", :verbose => false) }
    puts "SUCCESS. Failed when it was supposed to."
    puts ""
    
    puts "TEST: running '#{TESTGETOPTPROG}' with '--' to stop processing"
    assert_passes { sh("#{EXEC_PREFIX}#{TESTGETOPTPROG} -- monkey horse", :verbose => false) }
    puts "SUCCESS."
    puts ""
    
    puts "TEST: running '#{TESTGETOPTPROG}' with unknown option"
    assert_fails { sh("#{EXEC_PREFIX}#{TESTGETOPTPROG} -q", :verbose => false) }
    puts "SUCCESS. Failed when it was supposed to."
    puts ""
    
    puts "TEST: running '#{TESTGETOPTPROG}' with multiple known options"
    assert_passes { sh("#{EXEC_PREFIX}#{TESTGETOPTPROG} -b -f4 > opttestout", :verbose => false) }
    val = File.open("opttestout", "r") { |infile| infile.read }
    assert_passes { val =~ /You specified the -b option/ }
    assert_passes { val =~ /You specified the -f option with the parameter "4"/ }
    puts "SUCCESS."
    puts ""
  end
  
  desc "Run all tests"
  task :all => [:getopt, :ss, :ss2, :ss3] 
  
  desc "Run all sample_stats2 tests"
  task :ss2 => [:ss2vss, :ss2f]

  desc "Run all sample_stats3 tests"
  task :ss3 => [:ss3read, :ss3f]
end
