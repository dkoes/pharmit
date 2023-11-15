<?php
session_start();

require("lib.php");

if(isset($_REQUEST["logout"]))
{
	session_unset();
	session_destroy();
}

$gotinvalid = ""; //set to something if there was a botched login
if(isset($_REQUEST["login"]))
{
	$user = $_POST['user'];
	$pass = $_POST['pass'];

	$gotinvalid = login($user,$pass);
}

?>
<!DOCTYPE html>
<html>
	<head>
		<meta charset="UTF-8">
		<!--<meta name="viewport" content="width=960">!-->
		<title>pharmit: interactive exploration of chemical space</title>
		<link rel="stylesheet" type="text/css" href="index.css">
		<link href='http://fonts.googleapis.com/css?family=Open+Sans' rel='stylesheet' type='text/css'>
		<script src="js/jquery-2.1.3.js"  ></script>
		<script src="js/numeral.js" async></script>
		<script type="text/javascript" src="js/index.js" defer></script>
		<script type="text/javascript" src="js/msg.js" defer></script>
		<meta name="description" content="Interactive virtual screening of billions of structures using pharmacophore search molecular shapes.">
		<meta name="keywords" content="pharmacophores,pharmacophore search,virtual screening,molecular shape,energy minimization">

	</head>
	<body id="body">
		<div id="notice"> <div class=closebutton style="float:right"></div>
		<font color="red"> We are currently experiencing server issues (11/26). There is no ETA on a resolution. <font>
		</div>
		<div class="cont">
			<div class="section">
				<div class="cont-2">
					<div class="colwrapper">
						<div class="cont-3">
							<img src="images/pasted-image-448.png" alt="" class="img">
						</div>
					</div>
					<div class="colwrapper-2">
						<div class="cont-4">
							<div class="cont-5"><p class="para"><span class="font"><a href="https://github.com/dkoes/pharmit">code</a></span></p></div>
						</div>
						<div class="cont-6">
							<div class="cont-7"><p class="para-2"><span class="font-2"><a href="help.html">help</a></span></p></div>
						</div>
						<div class="cont-8">
							<div class="cont-9"><p class="para-3"><span class="font-3"><a href="mailto:dkoes@pitt.edu">contact</a></span></p></div>
						</div>
					</div>
				</div>
			</div>
		</div>
		<div class="section-2">
			<div class="cont-10">
				<div class="colwrapper-3">
					<div class="cont-11">
						<div class="cont-12"><p class="para-4"><span class="font-4">pharmit</span></p></div>
					</div>
					<div class="cont-13">
						<div class="cont-14"><p class="para-5"><span class="font-5">interactive exploration of chemical space</span></p></div>
					</div>
				</div>
			</div>
			<div class="cont-15"></div>
		</div>
		<div class="cont-16">
			<div class="section-3">
				<div class="cont-17">
					<div class="colwrapper-4">
						<div class="cont-18">
							<div class="cont-19"><p class="para-6"><span class="font-6">search</span></p></div>
						</div>
						<div class="cont-20">
							<div class="cont-21"><p class="para-7"><span class="font-7">virtual screening in your browser</span></p></div>
						</div>

						<div class="colwrapper-13">

							<div class="cont-27">
								<div class="cont-28"><p class="para-10"><span class="font-10"><a href="search.html">enter pharmit search</a></span></p></div>
							</div>
							<div class="cont-32">
								<div class="cont-33">
								<div class="pdbsearch">
								<form action='#' id="pdbform">
								start from PDB:
								<input type="text" id="pdbtext">
								<select disabled id="pdbligand">
								<option disabled="disabled" selected="true">ligand</option>
								</select>
								<br>
								binding site waters:
								<select id="pdbwater" >
								<option selected value="ignore">ignore</option>
								<option value="rec">treat as receptor</option>
								<option value="lig">treat as ligand</option>
								</select>
								<br>
								<input type="submit" class="pdbsubmit" id="pdbsubmit" value="submit" disabled>
								</form>
							</div>
								</div>
							</div>
							<div class="cont-45">
								<div class="cont-46"><p class="para-16"><span class="font-16"><a href="examples.html">examples</a></span></p></div>
							</div>
						</div>
					</div>
					<div class="colwrapper-5">
						<div class="cont-22">
							<div class="cont-23"><p class="para-8"><span class="font-8">create</span></p></div>
						</div>
						<div class="cont-24">
							<div class="cont-25"><p class="para-9"><span class="font-9">submit your own chemical libraries</span></p></div>
                            
<?php
if (!isset($_SESSION['userid']))
{
?>
			<div class="loginbox">
					<form action="index.php" method="POST">
					<input type="hidden" name="login" value="1">
					<div class="cont-30"><p class="para-11"><span class="font-11">log in to manage libraries</span></p></div>
					<div class="cont-36"><p class="para-13"><span class="font-13">email:</span></p></div>
					<input type="text" autofocus="autofocus" name="user" class="input-2" autocomplete="on">
					<div class="cont-40"><p class="para-14"><span class="font-14">password:</span></p></div>
					<input type="password" name="pass" class="input-2">
					<?php  echo("<div class=loginerror>$gotinvalid</div>"); ?>
					<p class="para-15"><span class="submit"><input type="submit" value="log in" class="submit" /></span></p></div>
					</form>
			</div>
			<div class="cont-48"><p class="para-17"><span class="font-17"><a href="create.php?op=register">register new account</a></span></p></div>
			<div class="cont-50"><p class="para-18"><span class="font-19"><a href="create.php?op=guestlogin">log in as guest</a></span></p></div>

<?php
}
else
{ //logged in
	$user = $_SESSION['userid'];
	$name = getname($user);
	$cnts = getcnts($user);
	$inprogress = $cnts["inprogress"];
	$completed = $cnts["completed"];
	echo("<div class=stuff>Welcome $name! You have $completed libraries available to search and $inprogress libraries under construction.</div>");
?>


			<div class="cont-48"><p class="para-17"><span class="font-17"><a href="index.php?logout=1">log out</a></span></p></div>
			<div class="cont-50"><p class="para-18"><span class="font-19"><a href="create.php?op=status">manage</a></span></p></div>
			<div class="cont-50"><p class="para-18"><span class="font-19"><a href="create.php">create</a></span></p></div>

<?php
}
?>
		</div>
		</div>
				</div>
			</div>
		</div>
		<div class="cont-51">
			<div class="section-4">
				<div class="cont-52">
					<div class="cont-53"><p class="para-19">
					<span class="font-21">pharmit currently has</span></p><p class="para-20">
					<span class="font-22"><span id="numberstandard"></span> built-in libraries </span>
					<span class="font-23"> containing </span>
					<span class="font-24"><span id="numconfsstandard"></span> conformations</span><span class="font-25"> of </span>
					<span class="font-26"><span id="nummolsstandard"></span> compounds</span><span class="font-27"> and</span></p>
					<p class="para-22"><span class="font-28"><span id="numberpublic"></span> publicly accessible user-contributed libraries </span>
					<span class="font-23"> containing</span></p><p class="para-21">
					<span class="font-24"><span id="numconfspublic"></span> conformations</span><span class="font-25"> of </span>
					<span class="font-26"><span id="nummolspublic"></span> compounds</span><span class="font-27">.</span></p>
					</p></div>
				</div>
			</div>
		</div>
<!--  TODO: implement dynamic updates
		<div class="cont-54"></div>
		<div class="section-5">
			<div class="cont-55">
				<div class="cont-56"><p class="para-23"><span class="font-29">updates</span></p><p class="para-24"><span class="font-30">Month Day, Year</span><span class="font-31"> PubChem rebuilt with xxxxxx total compounds</span></p><p class="para-25"><span class="font-32">Month Day, Year</span><span class="font-33"> Kinase Inhibitors created by Some Guy</span></p><p class="para-26"><span class="font-34">Month Day, Year</span><span class="font-35"> MolPort updated</span></p></div>
			</div>
		</div>
-->
		<div class="cont-57">
			<div class="section-6">
				<div class="cont-58">
					<div class="cont-59"><p class="para-27"><span class="font-36">acknowledgements</span></p>
					<p class="para-28"><span class="font-37">If you find this site useful please cite <span class="pharmitlink"><a href='http://nar.oxfordjournals.org/content/early/2016/04/19/nar.gkw287.long'>Pharmit: interactive exploration of chemical space</a><span></span>.</span></p></div>
					<p class="para-28"><span class="font-37">Pharmit is funded through R01GM108340 from the National Institute of General Medical Sciences and a generous gift from Relay Therapeutics. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institute of General Medical Sciences, the National Institutes of Health, or Relay Therapeutics.</span></p></div>
				</div>
			</div>
		</div>
		<div class="cont-60">
			<div class="section-7">
				<div class="cont-61">
					<div class="colwrapper-16">
						<div class="cont-62">
							<div class="cont-63"><p class="para-29"><span class="font-38">&COPY; David Ryan Koes, Charles Yuan and the University of Pittsburgh.</span></p></div>
						</div>
					</div>
					<div class="colwrapper-17">
						<div class="cont-64">
							<div class="cont-65"><p class="para-30"><span class="font-39"><a href="privacy.html">privacy policy</a></span></p></div>
						</div>
					</div>
				</div>
			</div>
		</div>
<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-73085320-1', 'auto');
  ga('send', 'pageview');

  $(document).ready(function() {
	  //$('#notice').slideDown();
	  $('#notice .closebutton').click(function() { $('#notice').slideUp(); });
  });
</script>		
	</body>
</html>
