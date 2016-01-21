<?php
session_start();

require("lib.php");

//produce html with error message
function failhtml($msg)
{
	echo("<html><title>Error</title><body>");
	printf("<h1>%s<h1>", $msg);
	echo("</body></html>");
	exit();
}

//from https://gist.github.com/tylerhall/521810
function generateStrongPassword($length = 9, $add_dashes = false, $available_sets = 'luds')
{
	$sets = array();
	if(strpos($available_sets, 'l') !== false)
		$sets[] = 'abcdefghjkmnpqrstuvwxyz';
	if(strpos($available_sets, 'u') !== false)
		$sets[] = 'ABCDEFGHJKMNPQRSTUVWXYZ';
	if(strpos($available_sets, 'd') !== false)
		$sets[] = '23456789';
	if(strpos($available_sets, 's') !== false)
		$sets[] = '!@#$%&*?';

	$all = '';
	$password = '';
	foreach($sets as $set)
	{
		$password .= $set[array_rand(str_split($set))];
		$all .= $set;
	}

	$all = str_split($all);
	for($i = 0; $i < $length - count($sets); $i++)
		$password .= $all[array_rand($all)];

	$password = str_shuffle($password);

	if(!$add_dashes)
		return $password;

	$dash_len = floor(sqrt($length));
	$dash_str = '';
	while(strlen($password) > $dash_len)
	{
		$dash_str .= substr($password, 0, $dash_len) . '-';
		$password = substr($password, $dash_len);
	}
	$dash_str .= $password;
	return $dash_str;
}

//html boilerplate for regular page
function headerhtml()
{
?>
<html>
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8">
<meta http-equiv="content-script-type" content="text/javascript">
<meta http-equiv="content-style-type" content="text/css">
<link rel="stylesheet" type="text/css" href="create.css" />
<script src="js/jquery-2.1.3.js"></script>
<link href='http://fonts.googleapis.com/css?family=Open+Sans' rel='stylesheet' type='text/css'>

<script>
function removeCheck(name, id) {
	//query user to make sure they want to remove database
	if(confirm("Are you sure you want to remove "+name+"?")) {
		location.replace("create.php?op=remove&id="+id);
	}
}
</script>

<title>pharmit library creation</title>
</head>

<body>

<?php
}

function footerhtml()
{
	echo("</body></html>");
}


//if we aren't logged in already, show login screen
if (!isset($_REQUEST["op"]) && !isset($_SESSION['userid']))
{
	headerhtml();
?>
<div class="cont">
<div class="cont-2">
<div class="loginpage">
<span class="font-3">log in</span><br>
<span class="font">log in to build or manage public or private libraries</span><br><br>
<div class="loginbox">
<form action="create.php" method="POST">
<input type="hidden" name="op" value="login">
<span class="font-2">email:</span>
<input type="text" autocomplete="on" autofocus="autofocus" name="user" size="60" class="input"><br><br>
<span class="font-2">password:</span>
<input type="password" name="pass" size="60" class="input"><br><br>
<input type="submit" value="log in" class="submit">
</form>
</div>
<br><br><br>
<span class="font-2">If you don't have an account or have lost your password you can <a href="create.php?op=register">register</a> for one
or <a href="create.php?op=guestlogin">log in as a guest</a>.</span>
<br><br><span class="font-2"><a href="index.php">return to pharmit</a></span><br>

</div></div></div>
<?php
	footerhtml();
}
else if(isset($_REQUEST["op"])) //operation
{
	$op = $_REQUEST["op"];
	switch ($op) {
		case "login":
			//check user/pass
			$user = $_POST['user'];
			$pass = $_POST['pass'];

			$err = login($user,$pass);
			if($err)
			{
				failhtml($err);
			}
			else
			{
				header("location:create.php");
				exit();
			}
			break;
		case "guestlogin":
			$_SESSION['userid']  = "guest";
			$_SESSION['maxprivatedbs'] = 0;
			$_SESSION['maxprivateconfs'] = 0;
			$_SESSION['maxconfs'] = 10000;
			header("location:create.php"); //reload now that session is set with no op
			exit();
			break;
		case "register":

			headerhtml();
			?>
			<div class="cont">
			<div class="cont-2">
			<div class="loginpage">
			<span class="font-3">register</span><br>
			<span class="font">provide your information and we will email you a password</span><br><br>
			<div class="loginbox">
			<form action="create.php" method="POST">
			<input type="hidden" name="op" value="doregister">
			<span class="font-2">email:</span>
			<input type="text" autofocus="autofocus" name="email" size="60" class="input"><br><br>
			<span class="font-2">name:</span>
			<input type="text" name="name" size="60" class="input"><br><br>
			<span class="font-2">institution:</span>
			<input type="text" name="place" size="60" class="input"><br><br>
			<input type="submit" value="submit registration" class="submit">
			</form>
			</div>
			<br><br><span class="font-2"><a href="index.php">return to pharmit</a></span><br>


			</div></div></div>
			<?php
				footerhtml();
			break;
		case "doregister":
			//create a user
			$email = $_POST['email'];
			$name = $_POST['name'];
			$place = $_POST['place'];

			if (!filter_var($email, FILTER_VALIDATE_EMAIL)) {
				failhtml("Invalid email: $email");
			}
			//generate a password - we assign a password and don't let the
			//user change since this way we (well, I) don't feel bad about
			//storing everything in clear text
			$pass = generateStrongPassword();

			//insert user (or replace, this is what you do when you lose your password)
			$db = new mysqli($db_host, $db_user, "", $db_name);
			if (mysqli_connect_errno())
				fail('MySQL connect', mysqli_connect_error());

			$stmt = $db->prepare('REPLACE INTO users (email, password, name, institution) VALUES (?,?,?,?)');
			$stmt->bind_param('ssss', $email,$pass,$name,$place) || fail('Bind register', $db->error);
			if (!$stmt->execute()) {
					failhtml('Unexpected error registering. Please try again later or contact the site administrator. ');
			}
			else {
				mail( $email, "Pharmit Password" ,
				"Your password is:\n$pass\n\nIf you lose this password you can simply re-register with the same email.",
						"From: do-not-reply@pharmit.csb.pitt.edu");
				echo("Your password has been mailed to $email.  Please check your spam filters. ");
				echo("<a href='create.php'>Continue</a>");
			}
			$stmt->close();

			break;
		case "remove":
			//get all the info for this user's databases
			$db = new mysqli($db_host, $db_user, "", $db_name);
			if (mysqli_connect_errno())
				fail('MySQL connect', mysqli_connect_error());
			
			($stmt = $db->prepare('UPDATE `databases` SET status="REMOVE" WHERE email=? AND id=?')) ||
			fail('Prepare databases', $db->error);
			$stmt->bind_param('ss', $_SESSION["userid"],$_REQUEST["id"]);
			if (!$stmt->execute()) {
				failhtml('Unexpected error removing library. Does it really belong to you?');
			}
			
			headerhtml();
			echo('<div class="cont"><div class="cont-2"><span class="font-3">manage</span><br>
					<span class="font">library removed</span><br><br>');
			footerhtml();
			
			break;
		case "status":
			//get all the info for this user's databases
			$db = new mysqli($db_host, $db_user, "", $db_name);
			if (mysqli_connect_errno())
				fail('MySQL connect', mysqli_connect_error());

			($stmt = $db->prepare('SELECT name, id, isprivate, status, message, submitted, completed, nummols, numconfs FROM `databases` WHERE email=? ORDER BY submitted DESC')) ||
				fail('Prepare databases', $db->error);
			$stmt->bind_param('s', $_SESSION["userid"]);
			if (!$stmt->execute()) {
				failhtml('Unexpected error checking status. ');
			}
			else {
				$stmt->store_result();
				headerhtml();
				echo('<div class="cont"><div class="cont-2"><span class="font-3">manage</span><br>
					<span class="font">view the status of your libraries</span><br><br>');

				if($stmt->num_rows > 0) { //have already created databases

					$stmt->bind_result($name, $id, $isprivate, $status, $message, $submitted, $completed, $nummols, $numconfs);
					while($stmt->fetch()) {
						if($status != "REMOVE") {
							echo('<div class="librarystatus"><span class="font-4">');
							echo("<b>$name</b>");
							echo(": $message <br>");
							if($isprivate) {
								echo("<b>Private</b><br>");
								echo("Access code: $id<br>");
							}
							else {
								echo("Public<br>");
							}
							echo("Submitted: $submitted <br>");
	
							if($status == "Completed") {
								echo("Completed: $completed<br>");
								echo(number_format($numconfs) . " conformers of ".number_format($nummols)." compounds");
								echo("<br><a class=\"removelink\" onclick=removeCheck(\"$name\",\"$id\")>Remove</a>");
							}
	
							echo("</div></span><br>");
						}
					}
				}
				else { //no databases
					echo('<span class="font-4">You have not created any databases.</span><br>');
				}
				echo('<br><span class="font-2"><a href="index.php">return to pharmit</a></span><br></div></div>');
				footerhtml();
			}


			break;
		case "logout":
			//remove session totally

			// Unset all of the session variables.
			session_unset();
			// If it's desired to kill the session, also delete the session cookie.
			// Note: This will destroy the session, and not just the session data!
			if (ini_get("session.use_cookies")) {
				$params = session_get_cookie_params();
				setcookie(session_name(), '', time() - 42000,
				$params["path"], $params["domain"],
				$params["secure"], $params["httponly"]
				);
			}

			// Finally, destroy the session.
			session_destroy();
			//back to login screen
			header("location:create.php");

			break;
	}
}
else //logged in, let's create some databases
{
	headerhtml();
	?>

	<div class="createpage">
	<div class="cont">
	<div class="cont-2">
	<span class="font-3">create</span><br>
	<div class="loginbox">
	<form id="createform" action="#" >
	<!--  TODO: add more input validation (name and file are required, check file name extension - sdf,smi,sdf.gz, or smi.gz - all in the client -->
	<input type="hidden" name="op" value="create">
	<span class="font">new databases from compounds</span><br><br>
	<span class="font-2">a short descriptive name of the database:</span><br>
	<input type="text" autofocus="autofocus" id="dbname" name="dbname" size="60" class="input-2"><br><br>
	<span class="font-2">a longer description: </span><br>
	<textarea rows="4" cols="50" id="description" name="description">
</textarea><br>
	<span class="font-4">please include any information you think may be useful, including contact information.</i></span><br><br>


  <?php
  //find out how many private databases this user already has
  $db = new mysqli($db_host, $db_user, "", $db_name);
  if (mysqli_connect_errno())
  	fail('MySQL connect', mysqli_connect_error());
  
  ($stmt = $db->prepare('SELECT COUNT(*) FROM `databases` WHERE email=? AND status != \'Error\' AND status != \'REMOVE\' AND isprivate != 0')) ||
  	fail('Prepare databases', $db->error);
  $stmt->bind_param('s', $_SESSION["userid"]);
  $stmt->execute();
  $stmt->store_result();
  $numprivate = 0;
  if($stmt->num_rows > 0) { //have valid username
  	$stmt->bind_result($numprivate) || fail('Bind cnt', $db->error);
  	if(!$stmt->fetch() && $db->errno)
  		fail('Fetch cnt', $db->error);
  }


  	echo('<span class="font-2">access:</span><br>');
  	if($numprivate > 0) echo("You have $numprivate private databases already built or pending.<br>");

 	echo('<select name="access">');
  	echo('   <option value="public">public - anyone will be able to view and search</option>');

 	$disabled = "";
	if($numprivate >= $_SESSION["maxprivatedbs"]) $disabled = " disabled ";
  	echo("<option $disabled value='private'>private - a passcode will be required to view and search</option>");
?>
	</select><br>
	<br><span class="font-2">compound file:</span>
	<input type="file" name="compounds">
	<br><br>
    <span class="font-4">Either .smi.gz or .sdf.gz.  Conformers will be automatically generated from SMILES files while the conformers
    present in the SDF file will be used.  Conformers of the same molecule are assumed to have the same name.  If the SMILES molecules
    do not have a name, they will be assigned an ID corresponding to the line number.</span><br>
    <br> <!-- want these messages to depend on the select above -->
<?php
  echo('<span class="font-4">');
  echo("You may create a maximum of ". number_format($_SESSION["maxprivatedbs"])." private databases each with at most ". number_format($_SESSION["maxprivateconfs"]) . " conformers.<br>");
  echo("Public databases may have as many as ".number_format($_SESSION["maxconfs"]) . " conformers. ");
  echo("These limits can be increased by submitting a short justification to dkoes@pitt.edu.<br>");
  echo("File sizes are limited to 200MB.  It is highly recommended that you submit a compressed (.gz) file.");
  echo("</span>");
  echo("<input type=hidden name=\"email\" value=\"".$_SESSION["userid"]."\">");

  ?>
  <br><br>
	<input type="submit"  id="submitbutton" value="submit" class="submit">
    </form>
   	<br><br>
   	<div id="createstatus" class="font-4"></div>
    <br><span class="font-2"><a href="index.php">return to pharmit</a></span><br>
	</div>


	</div></div></div>
	<script>
	var form = $('#createform').submit(function(event) {
		$('#createstatus').text("");

		if(!$('#dbname').val()) {
			console.log($('dbname').val());
			$('#createstatus').text("Missing database name.");
			return false;
		}
		if(!$('#description').val()) {
			$('#createstatus').text("Missing database description.");
			return false;
		}

		$('#createstatus').text("Uploading...");
		$('#submitbutton').prop('disabled',true);
		var fd = new FormData($('#createform').get(0));
		event.preventDefault(); //do our own submission with ajax
		$.ajax({
			url: '/fcgi-bin/createlib.fcgi',
			data: fd,
			cache: false,
			processData: false,
			contentType: false,
			type: 'POST'
		}).done(function(ret) {
			$('#submitbutton').prop('disabled',false);

			//this returns the unique id for referencing this library
			if(ret.lastIndexOf("Error",0) === 0) {
				$('#createstatus').html(ret);
			}
			else { //every okay so far, future errors will happen asynchronously
				$('#createstatus').html('Processing library.  Check status <a href="create.php?op=status">here</a>');
			}
		}).fail(function(x, status, e) {
			$('#createstatus').text("Error: "+e);
			$('#submitbutton').prop('disabled',false);

		});
		return false;
	});

	</script>
<?php

	footerhtml();
}
?>

