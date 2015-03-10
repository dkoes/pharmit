<?php 
session_start();

//configuration variable
$db_user = "pharmit";
#$db_host = "192.168.6.4";
$db_host = "localhost";
$db_name = "pharmit";
$debug = 1;


//this is called for errors that should not be exposed to the user
function fail($pub, $pvt = '')
{
	global $debug;
	$msg = $pub;
	if ($debug && $pvt !== '')
		$msg .= ": $pvt";

	//dkoes - try to only fail before generating html header..
	header('Content-Type: text/plain');

	exit("An error occurred ($msg).\n");
}

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
<script src="js/jquery-2.1.3.js" ></script> 

<title>Pharmit Library Creation</title>
</head>

<body>
			<form action="create.php" method="POST">
			<input type="hidden" name="op" value="logout">
			<input type="submit" value="Logout">
			</form>
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
<div class="loginpage">
Login to build or manage public or private libraries.
<div class="loginbox">
<form action="create.php" method="POST">
<input type="hidden" name="op" value="login">
Email:<br>
<input type="text" autocomplete="on" autofocus="autofocus" name="user" size="60"><br>
Password:<br>
<input type="password" autocomplete="on" name="pass" size="60"><br>
<input type="submit" value="Log in">
</form>
</div>

If you don't have an account or have lost your password you can <a href="create.php?op=register">register</a> for one
or you can <a href="create.php?op=guestlogin">login as guest</a>.

</div>
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

			$db = new mysqli($db_host, $db_user, "", $db_name);
			if (mysqli_connect_errno())
				fail('MySQL connect', mysqli_connect_error());
			
			($stmt = $db->prepare('SELECT password, maxprivatedbs, maxprivateconfs, maxdbs, maxconfs FROM users WHERE email=?')) ||
				fail('Prepare users', $db->error);
			$stmt->bind_param('s', $user) || fail('Bind user', $db->error);
			$stmt->execute();
			$stmt->store_result();
			if($stmt->num_rows > 0) { //have valid username
				$correctpass = "";
				$maxprivatedbs = 0; $maxprivateconfs = 0; $maxdbs = 0; $maxconfs = 0;
				$stmt->bind_result($correctpass, $maxprivatedbs, $maxprivateconfs, $maxdbs, $maxconfs) || fail('Bind pass', $db->error);
				if(!$stmt->fetch() && $db->errno)
					fail('Fetch pass', $db->error);
				
				if($correctpass == $pass) {
					session_regenerate_id();
					$_SESSION['userid']  = $user;
					$_SESSION['maxprivatedbs'] = $maxprivatedbs;
					$_SESSION['maxprivateconfs'] = $maxprivateconfs;
					$_SESSION['maxconfs'] = $maxconfs;
					session_write_close();
					header("location:create.php");
					exit();
			
				} else {
					failhtml("Invalid password for $user");
				}
			
			} else {
				failhtml("Invalid user");
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
			<div class="loginpage">
			Provide your email and information and we will email you a password.
			<div class="loginbox">
			<form action="create.php" method="POST">
			<input type="hidden" name="op" value="doregister">
			Email:<br>
			<input type="text" autofocus="autofocus" name="email" size="60"><br>
			Name:<br>
			<input type="text" name="name" size="60"><br>
			Institution:<br>
			<input type="text" name="place" size="60"><br>
			<input type="submit" value="Register">
			</form>
			</div>

			
			</div>
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
				echo("Your password has been mailed to $email.  Please check your spam filters.");
				echo("<a href='create.php'>Continue</a>");
			}
			$stmt->close();
			
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
				
				if($stmt->num_rows > 0) { //have already created databases
					
					$stmt->bind_result($name, $id, $isprivate, $status, $message, $submitted, $completed, $nummols, $numconfs);
					while($stmt->fetch()) {
						echo("<div class='librarystatus'>");
						echo("$name : $message <br>");
						echo("Submitted: $submitted <br>");
						if($isprivate) echo("<b>Private</b><br>");
						else echo("Public<br>");
						
						if($status == "Finished") {
							echo("Completed: $completed<br>");
							echo(number_format($numconfs) . " conformers of ".number_format($nummols));
						}
						
						echo("</div>");
					}
				} 
				else { //no databases
					echo("You have not created any databases.");
				}

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
	<div class="loginbox">
	<form id="createform" action="#" >
	
	<!--  TODO: add more input validation (name and file are required, check file name extension - sdf,smi,sdf.gz, or smi.gz - all in the client -->
	<input type="hidden" name="op" value="create">
	A short descriptive name of the database:<br>
	<input type="text" autofocus="autofocus" name="dbname" size="60"><br>
	A longer description. Please include any information you think may be useful.  May include html markup.	
	<textarea rows="4" cols="50" name="description">
	</textarea><br>
 <select name="access">
  <option value="public">Public - Anyone will be able to view and search</option>
  <option value="private">Private - A passcode is required to view and search.  There are additional limitations on the number and size of private databases. </option>
</select> 
<br> <!-- want these messages to depend on the select above -->
  <?php 
  echo("You may create a maximum of ". number_format($_SESSION["maxprivatedbs"])." private databases each with at most ". number_format($_SESSION["maxprivateconfs"]) . " conformers.");
  echo("Public databases may have as many as ".number_format($_SESSION["maxconfs"]) . " conformers. ");
  echo("These limits can be increased by submitting a short justification to dkoes@pitt.edu.<br>");
  echo("File sizes are limited to 200MB.  It is highly recommended that you submit a compressed (.gz) file.");
  
  echo("<input type=hidden name=\"email\" value=\"".$_SESSION["userid"]."\">");
  		
  ?>
	<br>Compound file.  Either .smi.gz or .sdf.gz.  Conformers will be automatically generated from SMILES files while the conformers
	present in the SDF file will be used.  Conformers of the same molecule are assumed to have the same name.  If the SMILES molecules
	do not have a name, they will be assigned an id corresponding to the line number.
	<input type="file" name="compounds">
	<br>
	<input type="submit" value="Submit">
	</form>
	<div id="createstatus"></div>
	</div>

	
	</div>
	<script>
	var form = $('#createform').submit(function(event) {
		$('#createstatus').text("Uploading...");
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
			//this returns the unique id for referencing this library
			if(ret.lastIndexOf("Error",0) === 0) {
				$('#createstatus').html(ret);
			}
			else { //every okay so far, future errors will happen asynchronously
				$('#createstatus').html('Processing library.  Check status <a href="create.php?op=status">here</a>');
			}
		}).fail(function(x, status, e) {
			$('#createstatus').text("Error: "+e);
			
		}); 
		return false;
	});

	</script>
<?php 
				
	footerhtml();
}
?>

